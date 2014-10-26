module Bio
  class Hamming
    NUCLEOTIDES = ['A', 'C', 'G', 'T']

    def distance(a, b)
      distance = 0
      for i in 0..a.length - 1 do
        distance += 1 if a[i] != b[i]
      end
      distance
    end

    def neighbors(pattern, d)
      return [pattern] if d == 0
      return NUCLEOTIDES if (pattern.length == 1)
      neighborhood = []
      suffix = pattern[1..pattern.length - 1]
      suffix_neighbors = neighbors(suffix, d)
      suffix_neighbors.each do |text|
        if (distance(suffix, text) < d)
          NUCLEOTIDES.each do |n|
            neighborhood.push(n +  text)
          end
        else
          neighborhood.push(pattern[0] + text)
        end
      end
      neighborhood
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
      frequencies.select.with_index do |p, i| 
        indexes.push(i) if p == count
      end
      indexes.each { |i| result.push(ori_c.number_to_pattern(i, k)) }
      result
    end


    def sorting_frequent_patterns(text, k, d)
      ori_c = OriC.new
      frequent_patterns = []
      neighborhoods = []
      index = []
      count = []

      for i in 0..(text.length - k + 1) do
        neighbors = neighbors(text[i..i + k], d)
        neighborhoods.push(neighbors)
      end

      interval = 0..(neighborhoods.length - 1)
      holder = neighborhoods.flatten
      for i in interval do
        pattern = holder[i]
        number =  ori_c.pattern_to_number(pattern)
        index[i] = number
        count[i] = 1
      end

      sorted = index.sort
      for i in interval do
        count[i + 1] = count[i] + 1 if (sorted[i] == sorted[i + 1])
      end

      for i in interval do
        if (count[i] == count.max)
          pattern = ori_c.number_to_pattern(sorted[i], k)
          frequent_patterns.push(pattern)
        end
      end

      frequent_patterns
    end
  end
end
