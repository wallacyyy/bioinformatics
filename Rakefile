require './lib/bio'

task :pattern_counter do
  sample = File.open('./lib/bio/chapter_1_part_2/data.txt', 'rb').read
  count = Bio::OriC.new(sample).count('AATCGTCAA')
  puts count   
end

task :most_frequents do
  sample = File.open('./lib/bio/chapter_1_part_2/data_2.txt', 'rb').read
  ori_c  = Bio::OriC.new(sample)
  puts ori_c.most_frequents(12).join(" ")
end
