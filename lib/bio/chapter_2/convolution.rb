module Bio
  class Convolution
    def peptide_table(m)
      "-#{m}"
    end

    def mass_table(p)
      p[1..-1].to_i
    end

    def spectral_convolution(spectrum)
      result = []
      for i in 0..spectrum.length - 1 do
        row = [spectrum[i]]
        for c in 0..i - 1
          mass = row[0] - spectrum[c]
          row.push(mass) if mass > 0
        end
        result.concat(row[1..row.length - 1])
      end
      result.sort
    end

    def restrict(spectrum, m)
      result = []
      c = spectral_convolution(spectrum)
      convolution = c.select { |p| p >= 57 and p < 200 }.sort.reverse
      return convolution if m >= convolution.length
      for i in m..convolution.length - 1 do
        if (convolution[i] < convolution[m - 1])
           return convolution[0..i - 1]
        end
      end
    end

    def convolution_cyclopeptide_sequence(spectrum, m, n)
      leaderboard = ['']
      leader_peptide = ''

      restricted = restrict(spectrum, m)

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
      leader_peptide[1..-1]
    end

    def cyclic_spectrum(peptide)
      spectrum = [0]
      prefix_mass = prefix_mass_value(peptide)
      peptide_mass = prefix_mass.last

      p = peptide.split('-')
      p.shift
      size = p.length

      for i in 0..size - 1 do
        for j in i + 1..size do
          spectrum.push(prefix_mass[j] - prefix_mass[i])
          if (i > 0 and j < size)
            spectrum.push(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
          end
        end
      end
      spectrum.sort
    end

    def linear_spectrum(peptide)
      spectrum = [0]
      prefix_mass = prefix_mass_value(peptide)

      p = peptide.split('-')
      p.shift
      size = p.length

      for i in 0..size - 1 do
        for j in i + 1..size do
          spectrum.push(prefix_mass[j] - prefix_mass[i])
        end
      end
      spectrum.sort
    end

    def prefix_mass_value(peptide)
      prefix_mass = [0]

      p = peptide.split('-').map(&:to_i)
      p.shift

      for i in 1..p.length do
        prefix_mass.push(prefix_mass[i - 1] + p[i - 1])
      end
      prefix_mass
    end

    def expand(spectrum, peptides)
      expanded = []

      if (peptides == [''])
        spectrum.each do |s|
          p = peptide_table(s)
          expanded.push(p) if p
        end
        return expanded
      end

      peptides.each do |p|
        for i in 57..200 do
          expanded.push(p + "-#{i}")
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
  end
end
