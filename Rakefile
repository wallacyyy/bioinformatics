require './lib/bio'

task :pattern_counter do
  sample = File.open('./lib/bio/chapter_1_part_2/data.txt', 'rb').read
  count = Bio::PatternCounter.new(sample, 'AATCGTCAA').count
  puts count   
end
