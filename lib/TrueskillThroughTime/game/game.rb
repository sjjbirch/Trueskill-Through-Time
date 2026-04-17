
module TrueskillThroughTime
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

    def graphical_model # Rubyists would crave a class for the graphical representation, ie the return of this.
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
        theta_div_pow2 = (sigma_trunc**2*d.sigma**2)/(d.sigma**2 - sigma_trunc**2) # the perfect crime
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

      # This usage of gr_tuple uses a fixed small epsilon because no expensive convergence message passing will occur.
      while gr_tuple(step, 1e-6) && (i < 10)
        step = [0.0, 0.0]
        (0...d.length-1).each do |e|
          d[e].prior = t[e].posterior_win - t[e+1].posterior_lose
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

  # Only used in Game
  class TeamVariable
    attr_accessor :prior, :likelihood_lose, :likelihood_win, :likelihood_draw
    def initialize(prior=NINF, likelihood_lose=NINF, likelihood_win=NINF, likelihood_draw=NINF)
      @prior = prior
      @likelihood_lose = likelihood_lose
      @likelihood_win = likelihood_win
      @likelihood_draw = likelihood_draw
    end

    # Renamed from p to prob. p is kernel print...
    # Never used, but a very useful method for building convenience methods.
    def prob
      @prior*@likelihood_lose*@likelihood_win*@likelihood_draw
    end

    def posterior_win
      @prior*@likelihood_lose*@likelihood_draw
    end

    def posterior_lose
      @prior*@likelihood_win*@likelihood_draw
    end

    # skill.likelihood called in batch, not of this class though
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
end
