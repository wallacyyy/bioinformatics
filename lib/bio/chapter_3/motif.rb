module Bio
  require 'matrix'

  class Motif
    ROWS = { 'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3 }

    def initialize
      @hamming = Hamming.new
      @origin = OriC.new
    end

    def enumeration(sample, k, d)
      k_mers = []
      enumerated = []
      patterns = []

      sample.each do |dna|
        k_mers.concat(dna.scan(/(?=(.{#{k}}))/))
      end

      k_mers.flatten.each do |p|
        patterns.concat(@hamming.neighbors(p, d))
      end

      patterns.each do |p|
        found = true
        sample.each do |dna|
          found = false unless @hamming.approximate_patterns(p, dna, d).any?
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
          d = @hamming.distance(pattern, p)
          hamming_distance = d if hamming_distance > d
        end
        distance += hamming_distance
      end
      distance
    end

    def median_string(sample, k)
      distance = Float::INFINITY
      for i in 0..4**k - 1 do
        pattern = @origin.number_to_pattern(i, k)
        d = distance_pattern_string(pattern, sample)
        if (distance >= d)
          distance = d
          median = pattern
        end
      end
      median
    end

    def matrix_profile(profile)
      matrix = Matrix.empty
      profile.each do |line|
       matrix = Matrix.rows(matrix.to_a << line.split.map(&:to_f)) 
      end
      matrix
    end

    def probability(pattern, profile)
      matrix = matrix_profile(profile)
      probability = 1
      for i in 0..pattern.length - 1
        p = matrix[ROWS[pattern[i]], i]
        probability *= p
      end
      probability.round(7)
    end

    def most_probable(dna, profile, k)
      patterns = dna.scan(/(?=(.{#{k}}))/).flatten
      max = 0
      result = ''
      patterns.each do |p|
        percent = probability(p, profile)
        if percent > max
          max = percent
          result = p
        end
      end
      result
    end
  end
end
