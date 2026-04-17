
module TrueskillThroughTime
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

    # Increase skill uncertainty based on passed dynamic uncertainty factor (usually global gamma)
    # occurs per match on non TTT histories, per unit time when times are provided based on elapsed time
    # Basically, when we do a rating event, jiggle the uncertainty to reflect reduced confidence in estimate from fading
    # and prevent settling into false minima/
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
end
