ENV['RACK_ENV'] = 'test'
require './lib/bio'
require 'benchmark'

RSpec.configure do |config|
  config.before { include Bio } 
  config.treat_symbols_as_metadata_keys_with_true_values = true
  config.run_all_when_everything_filtered = true
  config.filter_run :focus
  config.order = 'random'
end
