require 'spec_helper'

describe Bio::Motif do
  let(:motif) { Bio::Motif.new }

  context 'finds the motif enumeration' do
    it 'with a simple dataset' do
      sample = ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT']
      result = motif.enumeration(sample, 3, 1)
      expect(result.join(' ')).to eq('ATA ATT GTT TTT')
    end

    it 'with a complex dataset' do
      sample = []
      IO.foreach('./spec/fixtures/motif-enum-input.txt') do |line|
        sample.push(line.chomp)
      end

      result = motif.enumeration(sample, 5, 2)
      output = IO.read('./spec/fixtures/motif-enum-output.txt').chomp

      expect(result.sort).to eq(output.split.sort)
    end
  end

  it 'calculates the hamming distance between a pattern and a group of dna strings' do
    pattern = 'AAA'
    sample = 'TTACCTTAAC GATATCTGTC ACGGCGTTCG CCCTAAAGAG CGTCAGAGGT'.split
    expect(motif.distance_pattern_string(pattern, sample)).to eq(5)
  end


  context 'finds a median string' do
    it 'with a simple dataset' do
      dna = 'AAATTGACGCAT GACGACCACGTT CGTCAGCGCCTG GCTGAGCACCGG AGTACGGGACAG'.split
      expect(motif.median_string(dna, 3)).to eq('GAC')
    end
  end

  it 'finds the probability of a dna based on a profile' do
    pattern = 'TCGGGGATTTCC'
    profile = Bio::Reader.new('./spec/fixtures/profile.txt').read_lines
    expect(motif.probability(pattern, profile)).to eq(0.0205753)
  end

  context 'finds the most probable motif' do
    it 'with a simple dataset' do
      dna = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
      profile = Bio::Reader.new('./spec/fixtures/profile-2.txt').read_lines
      expect(motif.most_probable(dna, profile, 5)).to eq('CCGAG')
    end

    it 'with a complex dataset' do
      dna = IO.read('./spec/fixtures/profile-3-dna.txt').chomp
      profile = Bio::Reader.new('./spec/fixtures/profile-3.txt').read_lines
      expect(motif.most_probable(dna, profile, 6)).to eq('TGTCGC')
    end
  end

  it 'calculates a motif score' do
    profile = Bio::Reader.new('./spec/fixtures/motif-score.txt').read_lines
    expect(motif.score(profile)).to eq(30)
  end

  context 'search for a motif using greedy strategy' do
    it 'with a simple dataset' do
      dnas = Bio::Reader.new('./spec/fixtures/greedy-search.txt').read_lines
      result = motif.greedy_motif_search(dnas, 3, 5)
      expect(result).to eq(['CAG', 'CAG', 'CAA', 'CAA', 'CAA'])
    end

    it 'with a complex dataset' do
      dnas = Bio::Reader.new('./spec/fixtures/greedy-search-2-input.txt').read_lines
      output = Bio::Reader.new('./spec/fixtures/greedy-search-2-output.txt').read_lines
      result = motif.greedy_motif_search(dnas, 12, 25)
      expect(result).to eq(output)
    end
  end
end
