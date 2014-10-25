module Bio
  class OriC
    NUCLEOTIDES = ['A', 'T', 'C', 'G']

    def initialize(sample = '')
      @sample = sample
    end

    def nucleotides_table(k)
      @nucleotides || @nucleotides = NUCLEOTIDES.repeated_permutation(k)
                                     .map(&:join).sort
    end

    def count(pattern)
      @sample.scan(/(?=(#{pattern}))/).count
    end

    def most_frequents(k)
      frequencies = {}
      pattern = ""

      for i in 1..k do pattern << "." end
      candidates = @sample.scan(/(?=(#{pattern}))/).flatten.uniq
      candidates.each do |candidate| 
        frequencies["#{candidate}"] = count(candidate)
      end
      max = frequencies.values.max
      frequencies.select{ |k, v| v >= max }.keys.sort
    end

    def frequency_table(k, text = @sample)
      frequencies = []
      for i in 0..(4**k - 1) do frequencies[i] = 0 end
      for i in 0..(text.length - k) do
        max = i + k - 1
        pattern = text[i..max]
        j = pattern_to_number(pattern)
        frequencies[j] = frequencies[j] + 1
      end
      frequencies
    end   

    def clumps(k, l, t)
      code = { 'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3 }
      array_size = 4**k
      mask = array_size / 4
      frequent_patterns = Set.new
      frequencies = []
      index = 0
      queue = []

      for i in 0..(array_size - 1) do frequencies[i] = 0 end
      for j in 0..(k - 2) do
        index = (index % mask) * 4 + code[@sample[j]]
      end

      for j in (k - 1)..(@sample.length - 1) do
        index = (index % mask) * 4 + code[@sample[j]]
        queue.push(index) 

        if (j >= l)
          old_index = queue.shift
          frequencies[old_index] = frequencies[old_index] - 1
        end

        frequencies[index] = frequencies[index] + 1

        if (frequencies[index] >= t)
          pattern = number_to_pattern(index, k)
          frequent_patterns.add(pattern)
        end
      end
      frequent_patterns.to_a
    end

    def pattern_to_number(pattern)
      k = pattern.length
      nucleotides_table(k).index(pattern)
    end

    def number_to_pattern(index, k)
      table = nucleotides_table(k)
      table[index]
    end
  end
end
