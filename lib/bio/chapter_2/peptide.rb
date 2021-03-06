module Bio
  class Peptide
    PEPTIDE_TABLE =
      {  57 => ['G'], 71 => ['A'], 87 => ['S'], 97 => ['P'], 99 => ['V'],
         101 => ['T'], 103 => ['C'], 113 => ['I', 'L'], 114 => ['N'],
         115 => ['D'], 128 => ['K', 'Q'], 129 => ['E'], 131 => ['M'],
         137 => ['H'], 147 => ['F'], 156 => ['R'], 163 => ['Y'], 186 => ['W'] }

    MASS_TABLE =
      { 'G' => 57,  'A' => 71,  'S' => 87,  'P' => 97,  'V' => 99,
        'T' => 101, 'C' => 103, 'I' => 113, 'L' => 113, 'N' => 114,
        'D' => 115, 'K' => 128, 'Q' => 128, 'E' => 129, 'M' => 131,
        'H' => 137, 'F' => 147, 'R' => 156, 'Y' => 163, 'W' => 186 }

    def simple_score(peptide, spectrum, cyclic = true)
      cyclic ? f = cyclic_spectrum(peptide) : f = linear_spectrum(peptide)
      control = []
      score = 0
      f.each do |p|
        next if (control.count(p) > 0)
        a = f.count(p)
        b = spectrum.count(p)
        score += [a, b].min
        control.push(p)
      end
      score
    end

    def score(peptide, spectrum, cyclic = true)
      cyclic ? f = cyclic_spectrum(peptide) : f = linear_spectrum(peptide)
      score = 0
      count = 0
      length = spectrum.length - 1

      for i in 0..f.length - 1 do 
        next if (f[i] < spectrum[count])

        while (f[i] > spectrum[count] and count < length)
          count += 1
        end

        if (f[i] == spectrum[count] and count < length)
          score += 1
          count += 1
        end
      end
      score
    end

    def format_input(input)
      spectrum = input.split.map(&:to_i)
    end

    def is_consistent?(peptide, spectrum)
      linear = linear_spectrum(peptide)    
      consistent = true

      linear.each do |p|
        unless spectrum.include?(p)
          consistent = false
          break
        end
      end
      consistent
    end

    def format(peptide)
      daltons = []

      peptide.chars.to_a.each do |p|
        daltons.push(MASS_TABLE[p])
      end
      daltons.join('-') # \('-')/
    end

    def expand(spectrum, peptides)
      expanded = []

      if (peptides == [''])
        spectrum.each do |s|
          p = PEPTIDE_TABLE[s]
          expanded.push(p.first) if p
        end
        return expanded
      end

      peptides.each do |p|
        MASS_TABLE.keys.each do |m|
          next if (m == 'Q' or m == 'L')
          expanded.push(m + p)
        end
      end 
      expanded.uniq
    end

    def trim(leaderboard, spectrum, n)
      scores = {}

      for i in 0..leaderboard.length - 1 do
        peptide = leaderboard[i]
        scores[leaderboard[i]] = score(peptide, spectrum, false)
      end

      scores = scores.sort_by{|k,v| v}.reverse.to_h
      leader = scores.keys
      leader_scores = scores.values
      size = leader.length

      return leader if n >= size
      for i in n..leader.length - 1 do
        if (leader_scores[i] < leader_scores[n - 1])
          return leader[0..i - 1]
        end
      end
    end

    def leaderboard_sequence(spectrum, n)
      leaderboard = ['']
      leader_peptide = ''

      while leaderboard.any? do
        leaderboard = expand(spectrum, leaderboard)
        leaderboard.delete_if do |peptide|
          peptide_mass = prefix_mass_value(peptide).max
          spectrum_mass = spectrum.max
          if (peptide_mass == spectrum_mass)
            if (score(peptide, spectrum) > score(leader_peptide, spectrum))
              leader_peptide = peptide
              false
            end
          else 
            peptide_mass > spectrum_mass
          end
        end
        leaderboard = trim(leaderboard, spectrum, n)
      end
      format(leader_peptide)
    end

    def cyclopeptides_sequence(spectrum)
      peptides = ['']
      result = []

      while peptides.any? do
        peptides = expand(spectrum, peptides)
        peptides.delete_if do |peptide|
          mass = prefix_mass_value(peptide).max
          if (mass == spectrum.max)
            result.push(format(peptide)) if (cyclic_spectrum(peptide) == spectrum)
            true
          else
            !is_consistent?(peptide, spectrum)
          end
        end
      end
      result
    end

    def prefix_mass_value(peptide)
      prefix_mass = [0]

      for i in 1..peptide.length do
        prefix_mass.push(prefix_mass[i - 1] + MASS_TABLE[peptide[i - 1]])
      end
      prefix_mass
    end

    def cyclic_spectrum(peptide)
      spectrum = [0]
      prefix_mass = prefix_mass_value(peptide)
      peptide_mass = prefix_mass.last

      for i in 0..peptide.length - 1 do
        for j in i + 1..peptide.length do
          spectrum.push(prefix_mass[j] - prefix_mass[i])
          if (i > 0 and j < peptide.length)
            spectrum.push(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
          end
        end
      end
      spectrum.sort
    end

    def linear_spectrum(peptide)
      spectrum = [0]
      prefix_mass = prefix_mass_value(peptide)

      for i in 0..peptide.length - 1 do
        for j in i + 1..peptide.length do
          spectrum.push(prefix_mass[j] - prefix_mass[i])
        end
      end
      spectrum.sort
    end
  end
end
