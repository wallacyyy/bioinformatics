module Bio
  class Hamming
    NUCLEOTIDES = ['A', 'C', 'G', 'T']

    def frequent_sort(sample, k)
      ori_c = OriC.new
      frequents = []
      index = []
      count = []
      
      for i in 0..(sample.length - k)
        pattern = sample[i..i + k - 1]
        index[i] = ori_c.pattern_to_number(pattern)
        count[i] = 1
      end

      sorted_index = index.sort
      for i in 0..(sample.length - k)
        if (sorted_index[i] == sorted_index[i - 1])
          count[i] = count[i - 1] + 1
        end
      end

      max = count.max
      for i in 0..(sample.length - k)
        if (count[i] == max)
          pattern = ori_c.number_to_pattern(sorted_index[i], k)
          frequents.push(pattern)
        end
      end
      frequents
    end

    def frequent_patterns_sort(sample, k, d)
      ori_c = OriC.new
      index = []
      count = []
      frequents = []
      neighborhood = []

      for i in 0..(sample.length - k) do
        pattern = sample[i..i + k - 1]
        n = neighbors(pattern, d)
        neighborhood.concat(n)
      end

      neighborhood = neighborhood.sort
      interval = neighborhood.length - k
      for i in 0..interval
        pattern = neighborhood[i]
        index[i] = ori_c.pattern_to_number(pattern)
        count[i] = 1
      end

      sorted_index = index.sort
      for i in 0..interval
        if (sorted_index[i] == sorted_index[i - 1])
          count[i] = count[i - 1] + 1
        end
      end

      max = count.max
      for i in 0..interval
        if (count[i] == max)
          pattern = ori_c.number_to_pattern(sorted_index[i], k)
          frequents.push(pattern)
        end
      end
      frequents
    end

    def distance(a, b)
      distance = 0
      for i in 0..a.length - 1 do
        distance += 1 if a[i] != b[i]
      end
      distance
    end

    def distance_break(a, b, d)
      distance = 0
      for i in 0..a.length - 1 do
        distance += 1 if a[i] != b[i]
        break if (distance >= d)
      end
      distance
    end

    def neighbors(pattern, d)
      return [pattern] if (d == 0)
      return NUCLEOTIDES if (pattern.length == 1)
      neighborhood = []
      suffix = pattern[1..pattern.length - 1]
      suffix_neighbors = neighbors(suffix, d)
      suffix_neighbors.each do |text|
        if (distance_break(suffix, text, d) < d)
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
  end
end
