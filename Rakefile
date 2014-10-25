require './lib/bio'
require 'pry'

task :skew do
  sample = IO.read('./lib/bio/chapter_1/data/data_9.txt').chomp
  skew = Bio::Skew.new(sample)
  File.write('./result.txt', skew.minimal.join(' '))
end

task :pattern_counter do
  sample = ('./lib/bio/chapter_1/data/data.txt')
  count = Bio::OriC.new(sample).count('AATCGTCAA')
  puts count   
end

task :most_frequents do
  sample = IO.read('./lib/bio/chapter_1/data/data_2.txt')
  ori_c  = Bio::OriC.new(sample)
  puts ori_c.most_frequents(12).join(' ')
end

task :reverser do
  sample = IO.read('./lib/bio/chapter_1/data/ata_3.txt').chomp
  reversed = Bio::DnaComplement.new(sample).reverse
  File.write('./result.txt', reversed)
end

task :find_positions do
  sample = IO.read('./lib/bio/chapter_1/data/data_4.txt').chomp
  positions = Bio::DnaComplement.new(sample).positions('CAGTTCCCA')
  positions.each { |p| result << p.to_s + ' ' }
  File.write('./result.txt', result)
end

task :pattern_to_number do
  puts Bio::OriC.new.pattern_to_number('CTTCGC')
end

task :number_to_pattern do
  puts Bio::OriC.new.number_to_pattern(5437, 8)
end

task :frequency do
  sample = IO.read('./lib/bio/chapter_1/data/data_6.txt').chomp
  result =  Bio::OriC.new(sample).frequency_table(7)
  IO.write('./result.txt', result)
end

task :clumps do
  sample = IO.read('./lib/bio/chapter_1/data/data_8.txt').chomp
  clumps = Bio::OriC.new(sample).clumps(10, 506, 16)
  File.write('./result.txt', clumps.join(' '))
end
