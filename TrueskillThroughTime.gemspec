# frozen_string_literal: true

require_relative "lib/TrueskillThroughTime/version"

Gem::Specification.new do |spec|
  spec.name = "TrueskillThroughTime"
  spec.version = TrueskillThroughTime::VERSION
  spec.authors = ["Solomon Birch"]
  spec.email = ["birch.jj@gmail.com"]

  spec.summary = "A Ruby implementation of Trueskill Through Time."
  spec.homepage = "https://github.com/sjjbirch/Trueskill-Through-Time"
  spec.required_ruby_version = ">= 3.2"

  spec.metadata["homepage_uri"] = spec.homepage
  spec.metadata["source_code_uri"] = "https://github.com/sjjbirch/Trueskill-Through-Time"
  spec.metadata["changelog_uri"] = "https://github.com/sjjbirch/Trueskill-Through-Time/blob/master/CHANGELOG.md"

  # Specify which files should be added to the gem when it is released.
  # The `git ls-files -z` loads the files in the RubyGem that have been added into git.
  spec.files = Dir.chdir(__dir__) do
    Dir['{lib}/**/*', 'LICENSE', 'README.md']
  end
  spec.bindir = "exe"
  spec.executables = spec.files.grep(%r{\Aexe/}) { |f| File.basename(f) }
  spec.require_paths = ["lib"]

  # Uncomment to register a new dependency of your gem
  # spec.add_dependency "example-gem", "~> 1.0"

  # For more information and examples about making a new gem, check out our
  # guide at: https://bundler.io/guides/creating_gem.html
end
