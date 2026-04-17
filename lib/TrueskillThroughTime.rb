# frozen_string_literal: true

require_relative "TrueskillThroughTime/version"

module TrueskillThroughTime
  class Error < StandardError; end
  # Notes:
  #   Ported from Python trueskillthroughtime, v1.1.0
  #     Tuples are replaced with arrays.
  #       Functions that receive tuples are looking for quacks on [0] and [1].
  #         They will not fail if [>1] exists. But behaviour is undefined.
  #     Methods named p() are replaced with prob(). p() is a kernel function in Ruby.
  #     __rmul__ methods are replaced with coerce
  #     truediv methods are replaced with / (fdiv vs div in ruby etc)
  require 'json'
  require 'date'
  require_relative './TrueskillThroughTime/gaussian/gaussian'
  require_relative './TrueskillThroughTime/player/player'
  require_relative './TrueskillThroughTime/messages/messages'
  require_relative './TrueskillThroughTime/batch/batch'
  require_relative './TrueskillThroughTime/game/game'
  require_relative './TrueskillThroughTime/history/history'

  # Configuration Constants - defaults in case nothing is provided
  # Default mean average skill.
  # This is the value every player looks at and marvels at. Ironically it is a completely arbitrary and meaningless number.
  MU = 0.0              # Original Trueskill default is 25.
  # Class skill width. The difference in skill where a better player should beat a worse player 80% of the time.
  BETA = 1.0            #  Original Trueskill default is SIGMA/2
  # Default standard deviation of uncertainty. 95.4% of players should fall into +/- 2*SIGMA of the mean.
  # Note that this value only means anything in the context of BETA and GAMMA.
  SIGMA = BETA * 6.0    # Original Trueskill default is MU/3
  # Additive dynamics factor, also called tau in Trueskill (name clashes with Gaussian tau).
  GAMMA = BETA * 0.03   #  Original Trueskill default is BETA/100
  # The empirically observed probability of a global draw. If observation is not available, then this should be optimised.
  P_DRAW = 0.0
  # The thresholds for aborting convergence in Trueskill Through Time.
  # If the max changes in skill estimates from iterations of message passing are < EPSILON then convergence is achieved.
  EPSILON = 1e-6
  # Else stop after this number of iterations.
  ITERATIONS = 60

  # Derived Math Constants
  SQRT2 = Math.sqrt(2)
  SQRT2PI = Math.sqrt(2 * Math::PI)
  INF = Float::INFINITY

  N01 = Gaussian.new(0.0, 1.0)
  N00 = Gaussian.new(0.0, 0.0)
  NINF = Gaussian.new(0.0, INF)
  Nms = Gaussian.new(MU, SIGMA) # TODO: Check for deletion.
  #                                     Unused and potential source of bugs (uses module constant default MU/SIGMA).

  def erfc(x) # complementary error function, ie 1 - erf(x)
    # """(http://bit.ly/zOLqbc)"""
    z = x.abs
    t = 1.0 / (1.0 + z / 2.0)
    a = -0.82215223 + t * 0.17087277; b = 1.48851587 + t * a
    c = -1.13520398 + t * b; d = 0.27886807 + t * c; e = -0.18628806 + t * d
    f = 0.09678418 + t * e; g = 0.37409196 + t * f; h = 1.00002368 + t * g
    r = t * Math.exp(-z * z - 1.26551223 + t * h)
    if x.positive?
      r
    else
      2.0 - r
    end
  end

  # Approximate the normal inverse error function at y.
  def erfcinv(y)
    raise(ArgumentError, "Argument must be non-negative numbers") if y.negative?
    return -INF if y > 2
    return INF if y.zero?

    y = 2 - y if y >= 1

    t = Math.sqrt(-2 * Math.log(y / 2.0))
    x = -0.70711 * ((2.30753 + t * 0.27061) / (1.0 + t * (0.99229 + t * 0.04481)) - t)

    # Newton's method applied to approximate the normal inverse error function.
    3.times do
      err = erfc(x) - y
      x += err / (1.12837916709551257 * Math.exp(-(x**2)) - x * err)
    end

    y < 1 ? x : -x
  end

  # Given mu, sigma of a gaussian, compute tau and pi
  # Tau, pi being an alternative parameterisation of a gaussian using the precision.
  def tau_pi(mu, sigma)
    # Pi aka: Precision, inverse of variance. Not the 3.1415 etc constant.
    # Tau aka: Precision adjusted mean, ie mu * pi
    raise(ArgumentError, "Sigma must be greater than 0.") if (sigma + 1e-5) < 0.0

    if sigma > 0.0
      pi_ = sigma ** -2
      tau_ = pi_ * mu
    else
      pi_ = INF
      tau_ = INF
    end
    [tau_, pi_]
  end

  # Given tau, pi of a gaussian, compute mu, sigma
  def mu_sigma(tau_, pi_)
    raise(ArgumentError, "Sigma must be greater than 0.") if (pi_ + 1e-5) < 0.0

    if pi_ > 0.0
      sigma = Math.sqrt(1.0/pi_)
      mu = tau_ / pi_.to_f
    else
      sigma = INF
      mu = 0.0
    end
    [mu, sigma]
  end

  # cumulative density function
  def cdf(x, mu=0, sigma=1)
    0.5 * erfc(-(x - mu) / (sigma * SQRT2))
  end

  # Probability density function
  def pdf(x, mu, sigma)
    normaliser = (SQRT2PI * sigma)**-1
    functional = Math.exp(-((x - mu)**2) / (2*sigma**2))
    normaliser * functional
  end

  # percent point function
  def ppf(p, mu, sigma)
    mu - sigma * SQRT2 * erfcinv(2 * p)
  end

  # Basically, how much should we update values based on a win or loss with players of greater or lesser strength.
  # More abstractly, measures of "how much did we learn and in what direction" based on some outcome.
  def v_w(mu, sigma, margin, tie)
    # v aka: mean additive truncated gaussian function
    # w aka: variance multiplicative function
    if tie
      _alpha = (-margin-mu)/sigma
      _beta  = ( margin-mu)/sigma
      v = (pdf(_alpha, 0, 1)-pdf(_beta, 0, 1))/(cdf(_beta, 0, 1)-cdf(_alpha, 0, 1))
      u = (_alpha*pdf(_alpha, 0, 1)-_beta*pdf(_beta, 0, 1))/(cdf(_beta, 0, 1)-cdf(_alpha, 0, 1))
      w = - ( u - v**2 )
    else
      _alpha = (margin-mu)/sigma
      v = pdf(-_alpha, 0, 1) / cdf(-_alpha, 0, 1)
      w = v * (v + (-_alpha))
    end
    [v, w]
  end

  # Get the v and w, scale the incoming skill by them.
  # Used directly for a single game to directly return and update team posteriors. Otherwise used in approx.
  def trunc(mu, sigma, margin, tie)
    v, w = v_w(mu, sigma, margin, tie)
    mu_trunc = mu + sigma * v
    sigma_trunc = sigma * Math.sqrt(1-w)
    [mu_trunc, sigma_trunc]
  end

  # Call trunc, get back trunc'd mu and sigma based on v and w.
  # Return a gaussian.
  # Use the gaussian to update the win or loss likelihoods of a team.
  # Then to iterate that information to neighbouring nodes on the graph.
  def approx(n, margin, tie)
    mu, sigma = trunc(n.mu, n.sigma, margin, tie)
    Gaussian.new(mu, sigma)
  end

  def compute_margin(p_draw, sd)
    (ppf(0.5-p_draw/2.0, 0.0, sd)).abs
  end

  def podium(xs)
    sortperm(xs)
  end

  def sortperm(xs, reverse=false)
    sorted_indices = xs.each_with_index
                       .sort_by { |v, _| v }
                       .map { |_, i| i }
    sorted_indices.reverse! if reverse
    sorted_indices
  end

  # Tuplish methods
  def max_tuple(t1, t2) # TODO: Make class method of a tuple class )))))
    # t1,t2 are always Gaussians, t1 is always step
    # and this method is for making sure that the largest value is retained and propagates into step across iterations of a batch.
    # The step value is ultimately used to check for convergence against the Epsilon value by gr_tuple

    # Let's emulate python nan sorting behaviour!
    responses = [0.0, 0.0]
    responses[0] = if t1[0].to_f.nan? || t2[0].to_f.nan?
                     t1[0]
                   else
                     [t1[0], t2[0]].max
                   end
    responses[1] = if t1[1].to_f.nan? || t2[1].to_f.nan?
                     t1[1]
                   else
                     [t1[1], t2[1]].max
                   end
    responses
  end

  # Compare the values of a gaussian contained in a step variable to a threshold epsilon to see if
  # convergence has been achieved is accurate enough
  def gr_tuple(tup, threshold) # TODO: Make class method of a tuple class ((((
    tup.max > threshold
  end

  # Used to calculate the change in probability (ie gaussians) observed between runs of two batches
  def dict_diff(old, new)
    step = [0.0, 0.0]
    old.each_key do |a|
      step = max_tuple(step, old[a].delta(new[a])) # .delta being gaussian instance method
    end
    step
  end

  def get_composition(events)
    events.map {|event| event.teams.map {|team| team.items.map {|item| item.name} } }
  end

  def get_results(events) # TODO: Never called. Unused.
    events.map {|event| event.teams.map {|team| team.output } } # TODO: team.output usage 2.
  end

  def compute_elapsed(last_time, actual_time)
    if last_time == INF
      1
    elsif last_time == -INF
      0
    else
      actual_time - last_time
    end
  end

  def game_example(mu=MU, sigma=SIGMA, beta=BETA, p_draw=P_DRAW)
    a1 = Player.new(Gaussian.new(mu-500, sigma), beta=beta)
    a2 = Player.new(Gaussian.new(mu-500, sigma), beta=beta)
    a3 = Player.new(Gaussian.new(mu, sigma), beta=beta)
    a4 = Player.new(Gaussian.new(mu, sigma), beta=beta)
    a5 = Player.new(Gaussian.new(mu+5000, sigma), beta=beta)
    a6 = Player.new(Gaussian.new(mu+5000, sigma), beta=beta)
    team_a = [ a1, a2 ]
    team_b = [ a3, a4 ]
    team_c = [ a5, a6 ]
    teams = [team_a, team_b, team_c]
    result = [0.0, 1.0, 2.0]
    g = Game.new(teams, result, p_draw=p_draw)
    puts "Teams: #{g.teams}"
    puts "Likelihoods: #{g.likelihoods}"
    puts "Evidence: #{g.evidence}"
    puts "Posteriors: #{g.posteriors}"
    (0...teams.length).each do |i|
      puts "#{i}: #{g.performance(i)}"
    end
  end

  def history_example
    c1 = [["a"], ["b"]]
    c2 = [["b"], ["c"]]
    c3 = [["c"], ["a"]]
    composition = [c1, c2, c3]
    h = History.new(composition, gamma: 0.0)
    puts "#{h.learning_curves["a"]}"
    # [(1, N(mu=3.339, sigma=4.985)), (3, N(mu=-2.688, sigma=3.779))]
    puts "#{h.learning_curves["b"]}"
    # [(1, N(mu=-3.339, sigma=4.985)), (2, N(mu=0.059, sigma=4.218))]
    h.convergence
    puts "#{h.learning_curves["a"]}"
    # [[1, N(mu= 1.1776221018565704e-07, sigma= 2.3948083841791528)], [3, N(mu= -9.875327098012653e-08, sigma= 2.3948083506699978)]]
    puts "#{h.learning_curves["b"]}"
    # [[1, N(mu= -6.743217353120669e-08, sigma= 2.3948083803990317)], [2, N(mu= -6.888059735564812e-08, sigma= 2.394808375074641)]]
    puts "#{h.log_evidence}"
  end

  def liquid_example(path:"./example_data/liquid_t2_tournament_results.json")
    json_data = JSON.load(File.open(path))
    tournaments = {}
    json_data["tournament_positions"].each_pair do |k, v| # Symbolise keys and convert to dates, remove fake tournaments, remove tournaments with draws
      next unless v["results"].length > 2 # fake tournaments have < 3 teams
      results = {}
      had_ties = false
      v["results"].each_pair do |rk, rv| # numberise keys, then sort the results by them
        if rv.length > 1
          had_ties = true
          break
        else
          results[rk.to_f] = rv
        end
      end
      unless had_ties
        tournaments[k] = { :results => results.sort,
                           :day => DateTime.strptime(v["date"], "%Y-%m-%d %H:%M:%S").to_time.to_i / (60 * 60 * 24) }
      end
    end
    tournaments = tournaments.sort_by {|_, v| v[:day]}.to_h
    compositions = []
    results = []
    times = []
    tournaments.each_pair do |_, t_data|
      times << t_data[:day]
      t_results = []
      t_compositions = []
      t_data[:results].each_with_index do |(_, composition), i|
        t_results << t_data[:results].length - i
        t_compositions << composition.flatten
      end
      compositions << t_compositions
      results << t_results
    end

    h = History.new(compositions, results:results, times:times, mu:1200, sigma:400, beta:1320, gamma:19, p_draw: 0.01375)

    puts h.log_evidence
    h.convergence(epsilon: 0.01, iterations: 60, verbose: true)
    puts h.log_evidence

    h.evidence.each_with_index do |ev, idx|
      puts "#{idx}\t#{ev}\t#{tournaments.keys[idx]}"
    end

  end
end
