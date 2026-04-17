
module TrueskillThroughTime
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
end
