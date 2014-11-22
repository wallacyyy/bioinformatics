module Bio
  require 'matrix'
  require 'bigdecimal'
  require 'bigdecimal/util'

  class Motif
    ROWS = { 'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3 }

    def initialize
      @hamming = Hamming.new
      @origin = OriC.new
    end

    def init_frequencies 
      { 'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0 }
    end

    def greedy_motif_search(dnas, k, t)
      best_motifs = dnas.map{ |dna| dna[0, k] }
      first_dna_motifs = dnas.first.scan(/(?=(.{#{k}}))/).flatten

      first_dna_motifs.each do |motif|
        motifs = [motif]
        for i in 1..t - 1 do
          profile = to_profile(motifs) 
          most_probable_motif = most_probable(dnas[i], profile, k)
          #binding.pry if most_probable_motif == 'AGTGAATGTAAG'
          motifs.push(most_probable_motif)
        end
        best_motifs = motifs if score(motifs) <= score(best_motifs)
      end
      best_motifs
    end

    def to_profile(motifs)
      matrix = matrix_motif(motifs)
      profile = []

      matrix.column_vectors.each do |column|
        most_common = most_common_value(column) 
        frequencies = [0, 0, 0, 0]
        probabilities = []

        column.each do |p|
          frequencies[ROWS[p]] += 1 
        end

        total = frequencies.inject(0){ |t, element| t + element }

        frequencies.each_with_index do |v, i|
          p = (v.to_f/total).round(2)
          profile[i] ? profile[i].push(p) : profile[i] = [p]
        end
      end
      Matrix.rows(profile)
    end

    def most_common_value(array)
      array.group_by{ |e| e }.values.max_by(&:size).first
    end

    def score(motifs)
      matrix = matrix_motif(motifs)
      score = 0

      matrix.column_vectors.each do |column|
        most_common = most_common_value(column) 
        column.each do |p|
          score += 1 if p != most_common
        end
      end
      score
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

    def matrix_motif(motifs)
      return motifs if motifs.is_a?(Matrix)
      matrix = Matrix.empty

      motifs.each do |line|
       matrix = Matrix.rows(matrix.to_a << line.chars.to_a) 
      end
      matrix
    end

    def matrix_profile(profile)
      return profile if profile.is_a?(Matrix)
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
      probability
    end

    def most_probable(dna, profile, k)
      patterns = dna.scan(/(?=(.{#{k}}))/).flatten
      max = 0
      result = ''

      patterns.each do |p|
        percent = probability(p, profile)
        next if percent + max == 0 and not result.empty?
        if percent >= max
          max = percent
          result = p
        end
      end
      result
    end
  end
end
