# TrueskillThroughTime

This is a port of Trueskill Through Time (https://github.com/glandfried/TrueSkillThroughTime.py) to Ruby.

At present, v1.1.0 of TrueskillThroughTime.py is the ported version, corresponding to v1 of this gem.

In my test cases (large and complicated n vs n vs n etc FFA histories) it was approximately 1/4 to 1/3 faster than the target Python implementation, but I'm sure ymmv.

The gem was tested on Ruby 3.2.5 patch level 208, but there is no reason it shouldn't work on newer versions or older versions.

## Installation

Install the gem and add to the application's Gemfile by executing:

    $ bundle add Trueskill-Through-Time

## Usage

Usage patterns have been kept identical to those of the source version, so their documentation should basically work. Their documentation is at: https://trueskillthroughtime.readthedocs.io/en/latest/

If you do not converge your results in a History instance, then your usage is equivalent to the original Trueskill implementation, license and patents. If you do converge results, then your usage will fall under the later Trueskill Through Time implementation, license and patent.

If you want to go poking around, any time a Tuple is used in the Python version, a two-element Array is used in Ruby.

Finally, example usage is contained as the methods `game_example` and `history_example`, which mirror the examples given in the Python documentation. An additional, more complicated example is given as `liquid_example` using data from the Liquipedia API (https://liquipedia.net) under the CC-BY-SA 3 license. If you use this example data, then your usage must also comply with that license.

Otherwise, a convenience method `evidence` is added to the history class to get the discrete evidence by event. I will add some more convenience methods when I update the gem and alter the interfaces.

## Licensing
The original Trueskill and Trueskill Through Time are covered by the BSD license and separate Microsoft patents. One of these sets will apply to your usage, depending whether you converge your results. Notably, the only allowed uses either way are non-commercial products or XBox Live games.

The Liquipedia data used in the example is provided under CC-BY-SA 3.0, and if you use it then it will apply to your use.

Finally, this implementation is distributed under the Clear BSD License.

## Contributing

Bug reports and pull requests are welcome on GitHub at https://github.com/sjjbirch/TrueskillThroughTime. Probably hold off on PRs until it gets to v2, which I intend to be a less directly ported version with more convenience methods.
