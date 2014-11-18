require './lib/bio'
require 'pry'


task :enum do
  motif = Bio::Motif.new
  sample = []
  IO.foreach('./lib/bio/chapter_3/data/data.txt') do |line|
    sample.push(line.chomp)
  end
  result = motif.enumeration(sample, 5, 1)
  File.write('./result.txt', result.join(' '))
end


task :consistent do
  peptide = Bio::Peptide.new
  p = 'TCQ'
  spectrum = '0 71 99 101 103 128 129 199 200 204 227 230 231 298 303 328 330 332 333'
  puts peptide.is_consistent?(p, peptide.format_input(spectrum))
end

task :leader_convolution do
  peptide = Bio::Peptide.new
  convolution = Bio::Convolution.new
  spectrum = IO.read('./lib/bio/chapter_2/data/data_10.txt').chomp
  formatted = peptide.format_input(spectrum).sort
  m = 17
  n = 345
  result = convolution.convolution_cyclopeptide_sequence(formatted, m, n)
  puts result
end

task :convolution do
  peptide = Bio::Peptide.new
  c = Bio::Convolution.new
  spectrum = '0 57 118 179 236 240 301'
  formatted = peptide.format_input(spectrum) 
  result = c.spectral_convolution(formatted.sort)
  puts result
end

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
  sample = 'PEEP'
  spectrum = '0 97 97 129 129 194 203 226 226 258 323 323 323 355 403 452'
  peptide = Bio::Peptide.new
  formatted = peptide.format_input(spectrum) 
  puts peptide.simple_score(sample, formatted, false)
end

task :cyclo_sequence do
  input = IO.read('./lib/bio/chapter_2/data/data_3.txt').chomp
  peptide = Bio::Peptide.new
  result = peptide.cyclopeptides_sequence(peptide.format_input(input))
  File.write('./result.txt', result.join(' '))
end

task :spectrum do
  result = Bio::Peptide.new.cyclic_spectrum('TAIM')
  output = '0 71 101 113 131 184 202 214 232 285 303 315 345 416'
  f = result.join(' ')
  puts f == output
  puts f
  puts output
end

task :encode do
  sample = IO.read('./lib/bio/chapter_2/data/data_2.txt').chomp
  result = Bio::Rna.new.encode(sample, 'NHWHMLIM')
  File.write('./result.txt', result.join(' '))
end

task :rna_translate do
  sample = 'CCAAGAACAGAUAUCAAU'
  result = Bio::Rna.new.translate(sample)
  puts result
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
