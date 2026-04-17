
module TrueskillThroughTime

  class History
    attr_accessor :batches
    # Should be a class method of Agent,
    #     No, it should be ducktyped and live in History where it's called.
    def self.clean!(agents, last_time=false)
      agents.each_pair do |_, agent|
        agent.message = NINF
        agent.last_time = -INF if last_time
      end
    end

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
            a, # a unique player/agent/item/whatever identifier like a string player name taken from composition
            if priors.key?(a) # if we have a prior for the player then use that to instantiate an Agent
              Agent.new(priors[a], NINF, -INF)
            else
              Agent.new(Player.new(Gaussian.new(mu, sigma), beta, gamma), NINF, -INF) # if we don't then create a new one
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
      History.clean!(@agents)

      (0...@batches.length-1).to_a.reverse.each do |j|
        @batches[j+1].skills.each_key do |a|
          @agents[a].message = @batches[j+1].backward_prior_out(a)
        end
        old = @batches[j].posteriors # python copies
        # Recurses into iteration after populating new backward messages onto skills.backward from agents.message
        @batches[j].new_backward_info
        step = max_tuple(step, dict_diff(old, @batches[j].posteriors))
      end

      History.clean!(@agents)

      (1...@batches.length).each do |j|
        @batches[j-1].skills.each_key do |a|
          @agents[a].message = @batches[j-1].forward_prior_out(a)
        end
        old = @batches[j].posteriors # again, original copies
        # Recurses into iteration after populating new forward messages onto skills.forward from agents.receive
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

    def convergence(epsilon:EPSILON, iterations:ITERATIONS, verbose: false)
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

end
