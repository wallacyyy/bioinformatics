require './lib/bio'

task :pattern_counter do
  sample = File.open('./lib/bio/chapter_1/data.txt', 'rb').read
  count = Bio::OriC.new(sample).count('AATCGTCAA')
  puts count   
end

task :most_frequents do
  sample = File.open('./lib/bio/chapter_1/data_2.txt', 'rb').read
  ori_c  = Bio::OriC.new(sample)
  puts ori_c.most_frequents(12).join(" ")
end

task :reverser do
  sample = IO.readlines('./lib/bio/chapter_1/data_3.txt').first
  reversed = Bio::DnaComplement.new(sample).reverse
  File.write('./result.txt', reversed)
end

task :find_positions do
  sample = IO.readlines('./lib/bio/chapter_1/data_4.txt').first
  positions = Bio::DnaComplement.new(sample).positions('CAGTTCCCA')
  positions.each { |p| result << p.to_s + " " }
  File.write('./result.txt', result)
end

task :find_clumps do
  sample = IO.readlines('./lib/bio/chapter_1/data_5.txt').first
  clumps = Bio::DnaComplement.new(sample).clumps(10, 585, 20)
  result = ""
  clumps.each { |c| result << c.to_s + " " }
  File.write('./result.txt', result)
end
