module Bio
  class Skew
    def initialize(sample)
      @sample = sample
    end

    SUM_SKEW = { 'C' => -1, 
                 'G' => 1,
                 'A' => 0,
                 'T' => 0 }
    def find
      skew = [0]
      interval = 1..@sample.length
      for i in interval do
        skew[i] = skew[i - 1] + SUM_SKEW[@sample[i - 1]]
      end
      skew
    end

    def minimal
      skew = find
      min = skew.min
      result = []
      skew.each_index.select { |i| skew[i] == min }
    end
  end
end
