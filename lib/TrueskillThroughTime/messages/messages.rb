
module TrueskillThroughTime
  # Never instantiated or used.
  class DrawMessages # TODO: Function isn't clear. Never actually used in the original library either.
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

  # Only instantiated and used in Game
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
end
