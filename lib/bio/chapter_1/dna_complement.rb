module Bio
  class DnaComplement 
    attr_accessor :sample


    ASSOCIATIONS = { 'A' => 'T', 
                     'T' => 'A', 
                     'C' => 'G', 
                     'G' => 'C' } 

    def initialize(sample = '')
      @sample = sample
    end

    def reverse(sample = @sample)
      reversed = sample.reverse.gsub("\n", "")
      for i in 0..reversed.length - 1 do
        reversed[i] = ASSOCIATIONS[reversed[i]]
      end 
      reversed
    end

    def to_rna(sample = @sample)
      sample.gsub(/T/, 'U')
    end

    def positions(pattern, sample = @sample)
      sample.enum_for(:scan, /(?=(#{pattern}))/).map { Regexp.last_match.begin(0) }
    end
  end
end
