module Bio
  class Motif
    def hamming
      @hamming || @hamming = Hamming.new
    end

    def origin
      @origin || @origin = OriC.new
    end

    def enumeration(sample, k, d)
      k_mers = []
      enumerated = []
      patterns = []

      sample.each do |dna|
        k_mers.concat(dna.scan(/(?=(.{#{k}}))/))
      end

      k_mers.flatten.each do |p|
        patterns.concat(hamming.neighbors(p, d))
      end


      patterns.each do |p|
        found = true
        sample.each do |dna|
          found = false unless hamming.approximate_patterns(p, dna, d).any?
        end
        enumerated.push(p) if found
      end
      enumerated.uniq
    end

    def distance_pattern_string(pattern, sample)
      k = pattern.length
      distance = 0
      sample.each do |dna|
        hamming_distance = Float::INFINITY
        dna.scan(/(?=(.{#{k}}))/).flatten.each do |p|
          d = hamming.distance(pattern, p)
          hamming_distance = d if hamming_distance > d
        end
        distance += hamming_distance
      end
      distance
    end

    def median_string(sample, k)
      distance = Float::INFINITY
      for i in 0..4**k - 1 do
        pattern = origin.number_to_pattern(i, k)
        d = distance_pattern_string(pattern, sample)
        if (distance >= d)
          distance = d
          median = pattern
        end
      end
      median
    end
  end
end
