module TrueskillThroughTime

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
      # NOTE: This is the only used instantiation of Skill
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

    # TODO: Never called. Exists for the purpose of doing online rating, which I might add later.
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
        res[skill] = posterior(skill)
      end
      res
    end

    def within_prior(item)
      r = @agents[item.name].player
      r_p = (posterior(item.name) / item.likelihood)
      Player.new(Gaussian.new(r_p.mu, r_p.sigma), r.beta, r.gamma)
    end

    # Called only in iteration.
    def within_priors(event)
      @events[event].teams.map { |team| team.items.map { |item| within_prior(item) } }
    end

    def iteration(from=0)
      (from...@events.length).each do |e|
        teams = within_priors(e)
        result = @events[e].result
        weights = @events[e].weights
        # The advantage of instantiating this here is that it goes out of scope every loop,
        # the obvious disadvantage being that it's slower to create and tear down over and over.
        # We do become memory limited on large runs but the difference in size between a Game and an Event is not that big.
        # It might be worth folding all out Event usages into Games to reduce the number of classes.
        # Associated would be converting the Items to players
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

    # Skill instance message generator methods - could be Skill instance methods (would bloat their dispatch table tho)
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

  # Instantiated only in Batch at top level and in Batch in the unused add_events instance method
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

  class Agent # called only in History at top level
    attr_accessor :player, :message, :last_time
    def initialize(player, message, last_time)
      @player = player # a Player
      @message = message # a Gaussian, default NINF
      @last_time = last_time # a float, defaults to -infinity aka -INF
    end

    # Called once in Batch initialisation to create @skills Skill objects.
    # Called one other time in batch.new_forward_info instance method to update @skills Skill objects
    # Two additional usages in unused batch.add_events instance method
    def receive(elapsed)
      if @message != NINF # Real message.
        @message.forget(@player.gamma, elapsed)
      else # Instantiation/null state message
        @player.prior
      end
    end
  end

  # Called only in batch.
  # A lighterweight Agent used for marrying stats up to one later.
  # Always a component of a team.
  # Dumb name.
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

  # Instantiated only in Batch.
  # A holder for an array of teams and an array of results of those teams called output
  class Team
    attr_accessor :items, :output
    def initialize(items, output)
      @items = items # [Item, Item, Item]
      @output = output # [int, int, int] - inverse placements, ordinal list of "winningness" of that team by match
    end
  end

  # Instantiated only in Batch
  # Basically a lighter weight Game but in the context of a Batch, which is used solely for the purpose of Instantiating Games
  class Event
    attr_accessor :teams, :evidence, :weights
    def initialize(teams, evidence, weights)
      @teams = teams # [Team, Team, Team]
      @evidence = evidence # Float 0..1, likelihood of posterior given prior
      @weights = weights # [[float, float etc], [float, float etc], [float, float etc]] - partial play weights 0..1 for players on teams, can be []
      # length of weights == length of teams unless weights is empty
      # if weights is empty then the contribution of all players is assumed to be the full match
    end

    def to_s
      inspect
    end

    def inspect
      # { evidence: @evidence, teams: @teams, weights: @weights, result: result }.to_s
      "Event(#{names}, #{result})"
    end

    def names # Only called in inspect/to_s
      @teams.map { |team| team.items.map {|i| i.name} }
    end

    def result
      @teams.map { | team | team.output } # TODO: team.output usage 1.
    end
  end

end
