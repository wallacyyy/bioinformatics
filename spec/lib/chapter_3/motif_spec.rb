require 'spec_helper'

describe Bio::Motif do
  let(:motif) { Bio::Motif.new }

  xcontext 'finds the motif enumeration' do
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
end
