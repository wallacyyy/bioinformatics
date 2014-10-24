module Bio
  class OriC
    def initialize(sample)
      @sample = sample
    end

    def count(pattern)
      @sample.scan(/(?=(#{pattern}))/).count
    end

    def most_frequents(k)
      frequencies = {}
      pattern = ""

      for i in 1..k do pattern << "." end
      candidates = @sample.scan(/(?=(#{pattern}))/).flatten.uniq
      candidates.each { |candidate| frequencies["#{candidate}"] = count(candidate) }
      max = frequencies.values.max
      frequencies.select{ |k, v| v >= max }.keys.sort
    end
  end
end
