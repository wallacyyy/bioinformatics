module Bio
  class DnaComplement 
    ASSOCIATIONS = { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C' } 

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

    # How many strings with the length of k appears t times in an interval with length l
    def clumps(k, l, t)
      k_pattern = /(?=(.{#{k}}))/
      l_pattern = /(?=(.{#{l}}))/
      result = []

      intervals = @sample.gsub("\n", "").scan(l_pattern).flatten
      total = intervals.count
      intervals.each_with_index do |interval, index| 
        system('clear')
        puts "Reading interval #{index}/#{total}"
        found = interval.scan(k_pattern).flatten
        found.uniq.each do |p|
          matched = found.select { |e| e == p } 
          result.push(matched.first) if matched.count == t
        end
      end
      result.uniq
    end
  end
end
