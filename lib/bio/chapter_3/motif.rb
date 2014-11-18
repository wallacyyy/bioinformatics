module Bio
  class Motif
    def hamming
      @hamming || @hamming = Bio::Hamming.new
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
  end
end
