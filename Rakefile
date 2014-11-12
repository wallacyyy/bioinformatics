require './lib/bio'
require 'pry'

task :leaderboard do
  spectrum = IO.read('./lib/bio/chapter_2/data/data_6.txt').chomp
  peptide = Bio::Peptide.new
  formatted = peptide.format_input(spectrum) 
  puts peptide.leaderboard_sequence(formatted, 313)
end

task :trim do
  peptide = Bio::Peptide.new
  leaderboard = IO.read('./lib/bio/chapter_2/data/data_7.txt').chomp.split
  input = IO.read('./lib/bio/chapter_2/data/data_8.txt').chomp
  spectrum = peptide.format_input(input)
  result = peptide.trim(leaderboard, spectrum, 6) 
  File.write('./result.txt', result.join(' '))
end

task :score do
  sample = 'KAQSYPWQMEWPGQPEKKSIARWAWMDNSNTARNNIYPRH'
  spectrum = IO.read('./lib/bio/chapter_2/data/data_5.txt').chomp
  peptide = Bio::Peptide.new
  formatted = peptide.format_input(spectrum) 
  puts peptide.score(sample, formatted, false)
end

task :cyclo_sequence do
  input = IO.read('./lib/bio/chapter_2/data/data_3.txt').chomp
  peptide = Bio::Peptide.new
  result = peptide.cyclopeptides_sequence(peptide.format_input(input))
  File.write('./result.txt', result.join(' '))
end

task :spectrum do
  result = Bio::Peptide.new.cyclic_spectrum('PTYFRMWASRHFYA')
  File.write('./result.txt', result.join(' '))
end

task :encode do
  sample = IO.read('./lib/bio/chapter_2/data/data_2.txt').chomp
  result = Bio::Rna.new.encode(sample, 'NHWHMLIM')
  File.write('./result.txt', result.join(' '))
end

task :rna_translate do
  sample = IO.read('./lib/bio/chapter_2/data/data.txt').chomp
  result = Bio::Rna.new.translate(sample)
  File.write('./result.txt', result)
end

task :frequent_sort_reverse do
  sample = IO.read('./lib/bio/chapter_1/data/data_14.txt').chomp
  hamming = Bio::Hamming.new
  result = hamming.frequent_patterns_sort(sample, 10, 3, true)
  File.write('./result.txt', result.join(' '))
end

task :frequent_sort do
  sample = IO.read('./lib/bio/chapter_1/data/data_13.txt').chomp
  hamming = Bio::Hamming.new
  result = hamming.frequent_patterns_sort(sample, 10, 2)
  File.write('./result.txt', result.join(' '))
end

task :approximate do
  hamming = Bio::Hamming.new
  sample = 'CGTGACAGTGTATGGGCATCTTT'
  pattern = 'TGT'
  d = 1
  approximate = hamming.approximate_patterns(pattern, sample, d).count
  puts approximate
end

task :distance do
  data_1 = 'CTACAGCAATACGATCATATGCGGATCCGCAGTGGCCGGTAGACACACGT'
  data_2 = 'CTACCCCGCTGCTCAATGACCGGGACTAAAGAGGCGAAGATTATGGTGTG'
  hamming = Bio::Hamming.new
  distance = hamming.distance(data_1, data_2)
  puts distance
end

task :skew do
  sample = 'GATACACTTCCCAGTAGGTACTG'
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
  sample = 'CCAGATC'
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
  #clumps(k, l, t)
  sample = 'GCACAAGGCCGACAATAGGACGTAGCCTTGAAGACGACGTAGCGTGGTCGCATAAG' + 
           'TACAGTAGATAGTACCTCCCCCGCGCATCCTATTATTAAGTTAATT'
  clumps = Bio::OriC.new(sample).clumps(4, 30, 3)
  puts clumps
  File.write('./result.txt', clumps.join(' '))
end
