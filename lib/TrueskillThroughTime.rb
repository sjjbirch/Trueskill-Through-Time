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
  # Configuration Constants
  BETA = 1.0
  MU = 0.0
  SIGMA = BETA * 6
  GAMMA = BETA * 0.03
  P_DRAW = 0.0
  EPSILON = 1e-6
  ITERATIONS = 30

  # Derived Math Constants
  SQRT2 = Math.sqrt(2)
  SQRT2PI = Math.sqrt(2 * Math::PI)
  INF = Float::INFINITY

  def erfc(x)
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

  def erfcinv(y)
    raise(ArgumentError, "Argument must be non-negative numbers") if y.negative?
    return -INF if y > 2
    return INF if y.zero?

    y = 2 - y if y >= 1

    t = Math.sqrt(-2 * Math.log(y / 2.0))
    x = -0.70711 * ((2.30753 + t * 0.27061) / (1.0 + t * (0.99229 + t * 0.04481)) - t)

    3.times do
      err = erfc(x) - y
      x += err / (1.12837916709551257 * Math.exp(-(x**2)) - x * err)
    end

    y < 1 ? x : -x
  end

  def tau_pi(mu, sigma)
    raise(ArgumentError, "Sigma must be greater than 0.") if (sigma + 1e-5) < 0.0

    if sigma > 0.0
      pi_ = sigma ** -2
      tau_ = pi_ * mu
    else
      pi_ = INF
      tau_ = INF
    end
    [tau_, pi_] # python returns tuple
  end

  def mu_sigma(tau_, pi_)
    raise(ArgumentError, "Sigma must be greater than 0.") if (pi_ + 1e-5) < 0.0

    if pi_ > 0.0
      sigma = Math.sqrt(1.0/pi_)
      mu = tau_ / pi_.to_f
    else
      sigma = INF
      mu = 0.0
    end
    [mu, sigma] # again tuple in original
  end

  def cdf(x, mu=0, sigma=1)
    0.5 * erfc(-(x - mu) / (sigma * SQRT2))
  end

  def pdf(x, mu, sigma)
    normaliser = (SQRT2PI * sigma)**-1
    functional = Math.exp(-((x - mu)**2) / (2*sigma**2))
    normaliser * functional
  end

  def ppf(p, mu, sigma)
    mu - sigma * SQRT2 * erfcinv(2 * p)
  end

  def v_w(mu, sigma, margin, tie)
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
    [v, w] # again, a tuple
  end

  def trunc(mu, sigma, margin, tie)
    v, w = v_w(mu, sigma, margin, tie)
    mu_trunc = mu + sigma * v
    sigma_trunc = sigma * Math.sqrt(1-w)
    [mu_trunc, sigma_trunc]
  end

  def approx(n, margin, tie) # TODO: It seems like this should be a Guassian function
    mu, sigma = trunc(n.mu, n.sigma, margin, tie)
    # if $should_debug
    #   puts "approx(#{[n, margin, tie]}) (n, margin, tie) = mu #{mu} sigma #{sigma}"
    # end
    Gaussian.new(mu, sigma)
  end

  def compute_margin(p_draw, sd)
    (ppf(0.5-p_draw/2.0, 0.0, sd)).abs
  end

  def max_tuple(t1, t2)
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

  def gr_tuple(tup, threshold)
    tup.max > threshold
  end

  def podium(xs)
    sortperm(xs)
  end

  def sortperm(xs, reverse=false)
    sorted_indices = xs.each_with_index
                       .sort_by { |v, i| v }
                       .map { |v, i| i }
    sorted_indices.reverse! if reverse
    sorted_indices
  end

  def dict_diff(old, new)
    step = [0.0, 0.0]
    old.each_key do |a|
      step = max_tuple(step, old[a].delta(new[a])) # .delta is Guassian method....
    end
    step
  end

  class Gaussian
    extend Enumerable
    attr_accessor :mu, :sigma

    def initialize(mu=MU, sigma=SIGMA)
      raise(ArgumentError, "Sigma should be >= 0") unless sigma >= 0.0

      @mu = mu.to_f; @sigma = sigma.to_f
      @members = [@mu, @sigma]
    end

    def tau
      @sigma.positive? ? @mu * (@sigma**-2) : INF
    end

    def pi
      @sigma.positive? ? @sigma**-2 : INF
    end

    def each(&block)
      @members.each{|member| block.call(member)}
    end

    def to_s
      "N(mu= #{(@mu)}, sigma= #{(@sigma)})"
    end

    def inspect
      to_s
    end

    def +(other)
      unless other.is_a?(Gaussian)
        raise(ArgumentError, "Can only add gaussians to other gaussians - tried #{other.class}")
      end

      Gaussian.new(@mu + other.mu, Math.sqrt(@sigma**2 + other.sigma**2) )
    end

    def -(other)
      unless other.is_a?(Gaussian)
        raise(ArgumentError, "Can only subtract gaussians to other gaussians - tried #{other.class}")
      end

      Gaussian.new(@mu - other.mu, Math.sqrt(@sigma**2 + other.sigma**2) )
    end

    def *(other)
      if other.is_a?(Float)
        return NINF if other == INF

        # this is a gaussian instance, not -INF

        return Gaussian.new(other*@mu, other.abs*@sigma)

      elsif other.is_a?(Gaussian)
        if @sigma.zero? || other.sigma.zero?
          mu = other.mu / ((other.sigma**2/@sigma**2)+1)
          mu = @mu / ((@sigma**2 / other.sigma**2)+1) if @sigma.zero?
          sigma = 0.0
        else
          _tau = tau + other.tau
          _pi = pi + other.pi
          mu, sigma = mu_sigma(_tau, _pi)
        end
        return Gaussian.new(mu, sigma)
      else
        raise(ArgumentError, "Gaussian multiplication supports only floats or gaussian - tried #{other.class}")
      end
    end

    # Unused. Use / instead.
    def truediv(other)
      # Python equivalent of float division...
      # So for us just normal division:
    end

    def /(other)
      unless other.is_a?(Gaussian)
        raise(ArgumentError, "Gaussian division supports only gaussian - tried #{other.class}")
      end

      _tau = tau - other.tau
      _pi = pi - other.pi
      _mu, _sigma = mu_sigma(_tau, _pi)
      Gaussian.new(_mu, _sigma)
    end

    def rmul(other)
      # pythonism for other * this where mul is this * other
      # Ruby equivalent is to coerce. So:
    end

    def coerce(other)
      [self, other]
    end

    def forget(gamma, t)
      Gaussian.new(@mu, Math.sqrt(@sigma**2 + t*gamma**2))
    end

    def delta(other)
      [(@mu - other.mu).abs, (@sigma - other.sigma).abs]
    end

    def exclude(other)
      Gaussian.new(@mu - other.mu, Math.sqrt(@sigma**2 - other.sigma**2))
    end

    def isapprox?(other, tol=1e-4)
      ((@mu - other.mu).abs < tol) && ((@sigma - other.sigma).abs < tol)
    end
  end

  N01 = Gaussian.new(0.0, 1.0)
  N00 = Gaussian.new(0.0, 0.0)
  NINF = Gaussian.new(0.0, INF)
  Nms = Gaussian.new(MU, SIGMA)

  class Player
    attr_accessor :prior, :beta, :gamma, :prior_draw
    def initialize(prior=Gaussian.new(MU, SIGMA), beta=BETA, gamma=GAMMA, prior_draw=NINF)
      @prior = prior
      @beta = beta
      @gamma = gamma
      @prior_draw = prior_draw
    end

    def performance
      Gaussian.new(@prior.mu, Math.sqrt(@prior.sigma**2 + @beta**2))
    end

    def to_s
      "Player(Gaussian(mu=#{(@prior.mu)}, sigma=#{(@prior.sigma)}), beta=#{(@beta)}, gamma=#{(@gamma)})"
    end

    def inspect
      to_s
    end
  end

  class TeamVariable
    attr_accessor :prior, :likelihood_lose, :likelihood_win, :likelihood_draw
    def initialize(prior=NINF, likelihood_lose=NINF, likelihood_win=NINF, likelihood_draw=NINF)
      @prior = prior
      @likelihood_lose = likelihood_lose
      @likelihood_win = likelihood_win
      @likelihood_draw = likelihood_draw
    end

    # Renamed from p to prob. p is kernel print...
    def prob
      @prior*@likelihood_lose*@likelihood_win*@likelihood_draw
    end

    def posterior_win
      @prior*@likelihood_lose*@likelihood_draw
    end

    def posterior_lose
      @prior*@likelihood_win*@likelihood_draw
    end

    def likelihood
      @likelihood_win*@likelihood_lose*@likelihood_draw
    end

    def self.performance(team, weights) # I think this is used with weights as a single value
      weights = Array.new(team.length, weights) unless ((weights.is_a?(Array)) && (weights.length == team.length))
      res = N00
      team.zip(weights).each do |player, w|
        res += player.performance * w
      end
      res
    end

    def to_str
      to_s
    end

    def to_s
      inspect
    end

    def inspect
      "TTT::TeamVariable:#{self.object_id}:\t" + { prior: @prior, likelihood_lose: @likelihood_lose, likelihood_win: @likelihood_win, likelihood_draw: @likelihood_draw }.to_s
    end
  end

  class DrawMessages
    attr_accessor :prior, :prior_team, :likelihood_lose, :likelihood_win
    def initialize(prior=NINF, prior_team=NINF, likelihood_lose=NINF, likelihood_win=NINF)
      @prior = prior
      @prior_team = prior_team
      @likelihood_lose = likelihood_lose
      @likelihood_win = likelihood_win
    end

    # Renamed from p to prob. p() being a kernel method
    def prob
      @prior_team*@likelihood_lose*@likelihood_win
    end

    def posterior_win
      @prior_team*@likelihood_lose
    end

    def posterior_lose
      @prior_team*@likelihood_win
    end

    def likelihood
      @likelihood_win*@likelihood_lose
    end
  end

  class DiffMessages
    attr_accessor :prior, :likelihood
    def initialize(prior=NINF, likelihood=NINF)
      @prior = prior
      @likelihood = likelihood
    end

    def prob
      @prior*@likelihood
    end

    def to_str
      to_s
    end

    def to_s
      inspect
    end

    def inspect
      "TTT::DiffMessages:#{self.object_id}: \t" + {prior: @prior, likelihood: @likelihood}.to_s
    end
  end

  class Game
    attr_accessor :likelihoods, :teams, :evidence
    def initialize(teams, result=[], p_draw=0.0, weights=[])
      if (p_draw.negative? || p_draw > 1.0)
        raise(ArgumentError, "Draw probability #{p_draw} out of acceptable range: 1.0 > p_draw > 0.0.")
      end

      unless result.empty?
        if (result.length != teams.length)
          raise(ArgumentError, "Results provided but results for all teams not provided.")
        end
        if ((p_draw.zero?) && (result.length != result.uniq.length))
          raise(ArgumentError, "Draws provided in results but p_draw is 0.0.")
        end
      end
      unless weights.empty?
        if (weights.length != teams.length)
          raise(ArgumentError, "Weights provided but weights for all teams not provided.")
        end

        teams_without_weights = []
        teams.each_with_index do |team, idx|
          teams_without_weights << team if team.length != weights[idx].length
        end
        unless teams_without_weights.empty?
          raise(ArgumentError, "No weight provided for teams #{teams_without_weights}")
        end
      end

      @teams = teams
      @result = result # I hate that this isn't pluralised, but either is fine.
      @p_draw = p_draw
      weights = @teams.map {|team| Array.new(team.length, 1.0) } if weights.empty?
      @weights = weights
      @likelihoods = []
      @evidence = 0.0
      compute_likelihoods
    end

    # The number of teams in the game
    def length
      @teams.length
    end

    # The number of players in the Game
    def size
      i = 0 # should probably use inject and oneline it
      @teams.each do |team|
        i += team.length
      end
      i
    end

    def performance(i)
      # I'm not 100% how this is really meant to be used.
      # The semantics of Python zip() in python make it a little ambiguous.
      TeamVariable.performance(@teams[i], @weights[i])
    end

    def partial_evidence(d, margin, tie, e)
      mu = d[e].prior.mu
      sigma = d[e].prior.sigma
      @evidence *= if tie[e]
                     (cdf(margin[e], mu, sigma) - cdf(-margin[e], mu, sigma))
                   else
                     (1 - cdf(margin[e], mu, sigma))
                   end
    end

    def graphical_model
      if @result.empty?
        r = (0...@teams.length).to_a.reverse if @result.empty?
      else
        r = @result
      end

      o = sortperm(r, true)
      t = []
      d = []
      tie = []
      margin = []
      @teams.length.times do |i|
        t << TeamVariable.new(performance(o[i]), NINF, NINF, NINF)
      end

      (@teams.length-1).times do |i|
        d << DiffMessages.new(t[i].prior - t[i+1].prior, NINF)
        tie << (r[o[i]] == r[o[i+1]])
        if @p_draw.zero?
          margin << 0.0
        else
          sum = 0.0
          @teams[o[i]].each do |player|
            sum += player.beta**2
          end
          @teams[o[i+1]].each do |player|
            sum += player.beta**2
          end
          margin << compute_margin(@p_draw, Math.sqrt(sum))
        end
      end
      @evidence = 1.0
      [o, t, d, tie, margin]
    end

    def analytical_likelihood # "likelihood_analitico"
      o, t, d, tie, margin = graphical_model
      partial_evidence(d, margin, tie, 0)
      d = d[0].prior
      mu_trunc, sigma_trunc = trunc(d.mu, d.sigma, margin[0], tie[0])
      if d.sigma == sigma_trunc
        delta_div = d.sigma**2*mu_trunc - sigma_trunc**2*d.mu
        theta_div_pow2 = INF
      else
        delta_div = (d.sigma**2*mu_trunc - sigma_trunc**2*d.mu)/(d.sigma**2-sigma_trunc**2)
        delta_div_pow2 = (sigma_trunc**2*d.sigma**2)/(d.sigma**2 - sigma_trunc**2)
      end
      res = []
      t.length.times do |i|
        team = []
        @teams[o[i]].length.times do |j|
          mu == 0.0
          mu = @teams[o[i]][j].prior.mu + ( delta_div - d.mu)*(-1)**(i==1) unless d.sigma == sigma_trunc
          analytical_sigma = Math.sqrt(theta_div_pow2 + d.sigma**2 - @teams[o[i]][j].prior.sigma**2)
          team << Gaussian.new(mu, analytical_sigma)
        end
        res << team
      end

      if o[0] < o[1]
        [res[0], res[1]]
      else
        [res[1], res[0]]
      end
    end

    def likelihood_teams
      o, t, d, tie, margin = graphical_model

      step = [INF, INF]
      i = 0

      while gr_tuple(step, 1e-6) && (i < 10)
        step = [0.0, 0.0]
        (0...d.length-1).each do |e|
          d[e].prior = t[e].posterior_win - t[e+1].posterior_lose # prior giving wronnnggg
          partial_evidence(d, margin, tie, e) if i.zero?
          d[e].likelihood = approx(d[e].prior, margin[e], tie[e]) / d[e].prior
          likelihood_lose = t[e].posterior_win - d[e].likelihood
          step = max_tuple(step, t[e+1].likelihood_lose.delta(likelihood_lose))
          t[e+1].likelihood_lose = likelihood_lose
        end

        (1...d.length).to_a.reverse.each do |e|
          d[e].prior = t[e].posterior_win - t[e+1].posterior_lose
          partial_evidence(d, margin, tie, e) if ((i.zero?) && (e == (d.length-1)))
          d[e].likelihood = approx(d[e].prior, margin[e], tie[e]) / d[e].prior
          likelihood_win = t[e+1].posterior_lose + d[e].likelihood
          step = max_tuple(step, t[e].likelihood_win.delta(likelihood_win))
          t[e].likelihood_win = likelihood_win
        end
        i += 1
      end

      if d.length == 1
        partial_evidence(d, margin, tie, 0)
        d[0].prior = t[0].posterior_win - t[1].posterior_lose
        d[0].likelihood = approx(d[0].prior, margin[0], tie[0]) / d[0].prior
      end
      t[0].likelihood_win = t[1].posterior_lose + d[0].likelihood
      t[-1].likelihood_lose = t[-2].posterior_win - d[-1].likelihood
      (0...t.length).to_a.map { |e| t[o[e]].likelihood }
    end

    def compute_likelihoods
      weighted = 0
      @weights.each do |w|
        weighted += 1 if w != 1.0
      end
      if @teams.length > 2 || weighted.positive?
        m_t_ft = likelihood_teams
        res = []
        @teams.length.times do |e|
          e_likelihoods = []
          @teams[e].length.times do |i|
            weight_factor = if @weights[e][i] != 0.0
                              1 / @weights[e][i]
                            else
                              INF
                            end
            difference = m_t_ft[e] - performance(e).exclude(@teams[e][i].prior * @weights[e][i])
            e_likelihoods << (weight_factor * difference)
          end
          res << e_likelihoods
        end
        @likelihoods = res
      else
        @likelihoods = analytical_likelihood # untested
      end
    end

    def posteriors
      res = []
      @teams.each_with_index do |team, i|
        il = []
        team.each_with_index do |player, j|
          product = @likelihoods[i][j] * player.prior
          il << product
        end
        res << il
      end
      res
    end
  end

  class Skill
    attr_accessor :forward, :backward, :likelihood, :elapsed
    def initialize(forward=NINF, backward=NINF, likelihood=NINF, elapsed=0)
      @forward = forward
      @backward = backward
      @likelihood = likelihood
      @elapsed = elapsed
    end

    def to_s
      inspect
    end

    def inspect
      { 'forward'=>@forward, 'backward'=>@backward, 'likelihood'=>@likelihood, 'elapsed'=>@elapsed }.to_s
    end

  end

  class Agent
    attr_accessor :player, :message, :last_time
    def initialize(player, message, last_time)
      @player = player
      @message = message
      @last_time = last_time
    end

    def receive(elapsed)
      if @message != NINF
        @message.forget(@player.gamma, elapsed)
      else
        @player.prior
      end
    end
  end

  # Agents is a dict in python and their implementation is whack.
  def clean!(agents, last_time=false)
    agents.each_pair do |agent_label, agent|
      agent.message = NINF
      agent.last_time = -INF if last_time
    end
  end

  class Item
    attr_accessor :name, :likelihood
    def initialize(name, likelihood)
      @name = name
      @likelihood = likelihood
    end

    def to_s
      inspect
    end

    def inspect
      { 'name' => @name, 'likelihood' => @likelihood }.to_s
    end

  end

  class Team
    attr_accessor :items, :output
    def initialize(items, output)
      @items = items
      @output = output
    end
  end

  class Event
    attr_accessor :teams, :evidence, :weights
    def initialize(teams, evidence, weights)
      @teams = teams
      @evidence = evidence
      @weights = weights
    end

    def to_s
      inspect
    end

    def inspect
      # { evidence: @evidence, teams: @teams, weights: @weights, result: result }.to_s
      "Event(#{names}, #{result})"
    end

    def names
      @teams.map { |team| team.items.map {|i| i.name} }
    end

    def result
      @teams.map { | team | team.output }
    end
  end

  def get_composition(events)
    events.map {|event| event.teams.map {|team| team.items.map {|item| item.name} } }
  end

  def get_results(events)
    events.map {|event| event.teams.map {|team| team.output } }
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

  class Batch
    attr_accessor :skills, :events, :time
    def initialize(composition, results=[], time=0, agents={}, p_draw=0.0, weights=[])
      if (!results.empty? && (results.length != composition.length))
        raise(ArgumentError, "Results provided but composition.length != results.length")
      end
      if (!weights.empty? &&(weights.length != composition.length))
        raise(ArgumentError, "Weights provided but composition.length != weights.length")
      end

      this_agents = composition.flat_map { |teams| teams.flat_map { |team| team } }.to_set
      elapsed = this_agents.map { |a| [a, compute_elapsed(agents[a].last_time, time)] }.to_h
      @skills = this_agents.map { |a| [a, Skill.new(agents[a].receive(elapsed[a]), NINF, NINF, elapsed[a])] }.to_h
      @events = composition.each_with_index.map { |event_composition, e|
        Event.new(
          event_composition.each_with_index.map { |team_composition, t|
            Team.new(
              team_composition.each_with_index.map { |_, a|
                Item.new(composition[e][t][a], NINF)
              },
              results.empty? ? composition[e].length - t - 1 : results[e][t]
            )
          },
          0.0,
          weights.empty? ? weights : weights[e]
        )
      }
      @time = time
      @agents = agents
      @p_draw = p_draw
      iteration
    end

    def length
      @events.length
    end

    def add_events(composition, results=[])
      this_agents = composition.flat_map { |team| team.items }.to_set
      this_agents.each_with_index do |agent, i|
        elapsed = compute_elapsed(@agents[agent].last_time, @time)
        if @skills.include?(agent)
          @skills[agent] = Skill(@agents[agent].receive(elapsed), NINF, NINF, elapsed)
        else
          @skills[agent].elapsed = elapsed
          @skills[agent].forward = @agents[agent].receive(elapsed)
        end
      end

      _from = @events.length + 1

      composition.length.times do |e|
        @events << Event.new(
          composition[e].each_with_index.map { |team_composition, t|
            Team.new(
              team_composition.each_with_index.map { |_, a|
                Item.new(composition[e][t][a], NINF)
              },
              if results.empty?
                composition[e].length - t - 1
              else
                results[e][t]
              end
            )
          },
          0.0,
          if weights.nil? || weights.empty?
            weights
          else
            weights[e]
          end
        )
      end
      iteration(_from)
    end

    def posterior(agent)
      @skills[agent].likelihood * @skills[agent].backward * @skills[agent].forward
    end

    def posteriors
      res = {}
      @skills.each_key do |skill|
        res[skill] = posterior(skill) # I think this is fucked because posterior doesn't index this way
      end
      res
    end

    def within_prior(item)
      r = @agents[item.name].player
      r_p = (posterior(item.name) / item.likelihood)
      Player.new(Gaussian.new(r_p.mu, r_p.sigma), r.beta, r.gamma)
    end

    def within_priors(event)
      @events[event].teams.map { |team| team.items.map { |item| within_prior(item) } }
    end

    def iteration(from=0)
      (from...@events.length).each do |e|
        teams = within_priors(e)
        result = @events[e].result
        weights = @events[e].weights
        g = Game.new(teams, result, @p_draw, weights)
        @events[e].teams.each_with_index do |team, t|
          team.items.each_with_index do |item, i|
            @skills[item.name].likelihood = (@skills[item.name].likelihood / item.likelihood) * g.likelihoods[t][i]
            item.likelihood = g.likelihoods[t][i]
          end
        end
        @events[e].evidence = g.evidence
      end
    end

    def convergence(epsilon=1e-6, iterations = 20)
      step = [INF, INF]
      i = 0
      while gr_tuple(step, epsilon) && (i < iterations)
        old = posteriors # python makes a copy
        iteration
        step = dict_diff(old, posteriors)
        i += 1
      end
      i
    end

    def forward_prior_out(agent)
      @skills[agent].forward * @skills[agent].likelihood
    end

    def backward_prior_out(agent)
      n = @skills[agent].likelihood * @skills[agent].backward
      n.forget(@agents[agent].player.gamma, @skills[agent].elapsed)
    end

    def new_backward_info
      @skills.each_key do |a|
        @skills[a].backward = @agents[a].message
      end
      iteration # python explicitly: return self.iteration() -> and iteration() returns None???
      nil
    end

    def new_forward_info
      @skills.each_key do |a|
        @skills[a].forward = @agents[a].receive(@skills[a].elapsed)
      end
      iteration # python explicitly: return self.iteration() -> and iteration() returns None?
      nil
    end

    def to_s
      "Batch(time=#{@time}, events=#{@events})"
    end

    def inspect
      to_s
    end

  end

  class History
    attr_accessor :batches

    def initialize(composition, results:[], times:[], priors:{}, mu:MU, sigma:SIGMA, beta:BETA, gamma:GAMMA, p_draw:P_DRAW, weights:[])
      if results.length.positive? && results.length != composition.length
        raise(ArgumentError, "Results provided but composition.length (#{composition.length}) != results.length} (#{results.length})")
      end
      if times.length.positive? && times.length != composition.length
        raise(ArgumentError, "Times provided but composition.length (#{composition.length}) != times.length} (#{times.length})")
      end
      if weights.length.positive? && weights.length != composition.length
        raise(ArgumentError, "Weights provided but composition.length (#{composition.length}) != weights.length} (#{weights.length})")
      end

      @size = composition.length
      @batches = []
      @agents = Hash[
        composition.flatten(2).uniq.map do |a|
          [
            a,
            if priors.key?(a)
              Agent.new(priors[a], NINF, -INF)
            else
              Agent.new(Player.new(Gaussian.new(mu, sigma), beta, gamma), NINF, -INF)
            end
          ]
        end
      ]
      @mu = mu
      @sigma = sigma
      @p_draw = p_draw
      @time = times.length.positive?
      trueskill(composition, results, times, weights)

    end

    def trueskill(composition, results, times, weights)
      o = times.length.positive? ? sortperm(times) : (0...composition.length).to_a
      i = 0
      while i < @size
        j = i + 1
        t = times.length.zero? ? j : times[o[i]]
        while (times.length.positive?) && (j < @size) && (times[o[j]] == t) # I don't even
          j += 1
        end
        b = if results.length.positive?
              Batch.new(
                o[i...j].map { |k| composition[k] },
                o[i...j].map { |k| results[k] },
                t,
                @agents,
                @p_draw,
                weights.length.positive? ? o[i...j].map { |k| weights[k] } : weights
              )
            else
              Batch.new(
                o[i...j].map { |k| composition[k] },
                [],
                t,
                @agents,
                @p_draw,
                weights.length.positive? ? o[i...j].map { |k| weights[k] } : weights
              )
            end
        @batches << b
        b.skills.each_key do |a|
          @agents[a].last_time = (@time ? t : INF)
          @agents[a].message = b.forward_prior_out(a)
        end
        i = j
      end
    end

    def iteration
      step = [0.0, 0.0]
      clean!(@agents)

      (0...@batches.length-1).to_a.reverse.each do |j|
        @batches[j+1].skills.each_key do |a|
          @agents[a].message = @batches[j+1].backward_prior_out(a)
        end
        old = @batches[j].posteriors # python copies
        @batches[j].new_backward_info
        step = max_tuple(step, dict_diff(old, @batches[j].posteriors))
      end

      clean!(@agents)

      (1...@batches.length).each do |j|
        @batches[j-1].skills.each_key do |a|
          @agents[a].message = @batches[j-1].forward_prior_out(a)
        end
        old = @batches[j].posteriors # again, original copies
        @batches[j].new_forward_info
        step = max_tuple(step, dict_diff(old, @batches[j].posteriors))
      end

      if @batches.length == 1
        old = @batches[0].posteriors
        @batches[0].convergence
        step = max_tuple(step, dict_diff(old, @batches[0].posteriors))
      end

      step
    end

    def convergence(epsilon:EPSILON, iterations:ITERATIONS, verbose:true)
      step = [INF, INF]
      i = 0

      while gr_tuple(step, epsilon) && (i < iterations)
        step = iteration
        puts "Iteration = #{i} , #{step}" if verbose
        i += 1
      end
      puts "End" if verbose

      [step, i]
    end

    def learning_curves
      res = {}
      @batches.each do |batch|
        batch.skills.each_key do |a|
          res[a] ||= []
          res[a] << [batch.time, batch.posterior(a)]
        end
      end
      res
    end

    def evidence
      ev = []
      @batches.each do |batch|
        batch.events.each do |event|
          ev << event.evidence
        end
      end
      ev
    end

    def log_evidence
      # Math.log(@batches.flat_map { |b| b.events }.map { |event| event.evidence }.reduce(1, :*))
      @batches.flat_map { |b| b.events }.map { |event| Math.log(event.evidence) }.sum
    end

    def to_s
      "History(Events=#{@size}, Batches=#{@batches.length}, Agents=#{@agents})"
    end

    def inspect
      to_s
    end

    def length
      @size
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

  def history_example(mu=MU, sigma=SIGMA, beta=BETA, p_draw=P_DRAW)
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
    puts "Foo"
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
    tournaments = tournaments.sort_by {|k, v| v[:day]}.to_h
    # Need to get to compositions, results, times: results needs to be inverse of position
    compositions = []
    results = []
    times = []
    tournaments.each_pair do |tournament_name, t_data|
      times << t_data[:day]
      t_results = []
      t_compositions = []
      t_data[:results].each_with_index do |(position, composition), i|
        t_results << t_data[:results].length - i
        # t_compositions << composition.map {|player| [player]}
        t_compositions << composition.flatten
      end
      compositions << t_compositions
      results << t_results
    end

    h = History.new(compositions, results:results, times:times, mu:1200, sigma:400, beta:1320, gamma:19, p_draw: 0.01375)

    puts h.log_evidence
    h.convergence(epsilon: 0.01, iterations: 60, verbose: true)
    puts h.log_evidence
    puts h.learning_curves.first

    h.evidence.each_with_index do |ev, idx|
      puts "#{idx}\t#{ev}\t#{tournaments.keys[idx]}"
    end

  end
end
