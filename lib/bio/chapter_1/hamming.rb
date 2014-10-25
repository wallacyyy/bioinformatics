module Bio
  class Hamming
    def distance(a, b)
      distance = 0
      for i in 0..a.length - 1 do
        distance += 1 if a[i] != b[i]
      end
      distance
    end

    def approximate_patterns(pattern, sample, d)
      i = 0
      k = pattern.length
      result = []
      while i < (sample.length - k + 1) do
        interval = i..(i + k - 1)
        mask = sample[interval]
        result.push(i) if (distance(mask, pattern) <= d)
        i += 1
      end
      result
    end

    def frequent_patterns(sample, k, d)
      ori_c = OriC.new
      table = ori_c.nucleotides_table(k)
      frequencies = []
      indexes = []
      result = []

      for a in 0..(table.length - 1) do frequencies[a] = 0 end

      max = sample.length - k
      for i in 0..max do
        system('clear')
        puts "#{i}/#{max}"
        interval = i..(i + k - 1)
        mask = sample[interval]
        table.each_with_index do |pattern, index|
          frequencies[index] += 1 if (distance(mask, pattern) <= d)
        end
      end
      count = frequencies.max
      frequencies.select.with_index { |p, i| indexes.push(i) if p == count }
      indexes.each { |i| result.push(ori_c.number_to_pattern(i, k)) }
      result
    end
  end
end
