module Bio
  class DnaComplement 
    ASSOCIATIONS = { 'A' => 'T', 
                     'T' => 'A', 
                     'C' => 'G', 
                     'G' => 'C' } 

    def initialize(sample)
      @sample = sample
    end

    def reverse
      reversed = @sample.reverse.gsub("\n", "")
      for i in 0..reversed.length - 1 do
        reversed[i] = ASSOCIATIONS[reversed[i]]
      end 
      reversed
    end

    def positions(pattern)
      @sample.enum_for(:scan, /(?=(#{pattern}))/).map { Regexp.last_match.begin(0) }
    end
  end
end