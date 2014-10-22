module Bio
  class PatternCounter
    def initialize(sample, pattern)
      @sample = sample
      @pattern = pattern
    end

    def count
      @sample.scan(/(?=(#{@pattern}))/).count
    end
  end
end
