require 'spec_helper'

describe Bio::Hamming do
  let(:hamming) { Bio::Hamming.new }

  it 'calculates the hamming distance between two dna strings' do
    distance = hamming.distance('GGGCCGTTGGT', 'GGACCGTTGAC')
    expect(distance).to eq(3)
  end

  it 'finds the approximate patterns based on a hamming distance' do
    sample = 'CGCCCGAATCCAGAACGCATTCCCATATTTCG' + 
             'GGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT'
    pattern = 'ATTCTGGA'
    result = hamming.approximate_patterns(pattern, sample, 3) 
    expect(result.join(' ')).to eq('6 7 26 27')
  end

  it 'finds the number of approximate patterns based on a hamming distance' do
    sample = 'TTTAGAGCCTTCAGAGG'
    pattern = 'GAGG'
    result = hamming.approximate_patterns(pattern, sample, 2).count
    expect(result).to eq(4)
  end

  context 'finds most frequent patterns based on a hamming distance' do
    it 'runs with a simple dna string' do
      sample = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
      result = hamming.frequent_patterns(sample, 4, 1)
      expect(result.join(' ')).to eq('ATGC ATGT GATG')
    end

    xit 'runs with a complex dna string' do
      sample = IO.read('./spec/fixtures/mismatch.txt').chomp
      result = hamming.frequent_patterns(sample, 10, 2)
      expect(result.join(' ')).to eq('GCACACAGAC GCGCACACAC')
    end
  end

  it 'finds the neighbord patterns' do
    pattern = 'ACG'
    result = hamming.neighbors(pattern, 1)
    output = ['CCG', 'TCG', 'GCG', 'AAG', 'ATG', 'AGG', 
              'ACA', 'ACC', 'ACT', 'ACG'].sort
    expect(result.sort).to eq(output)
  end
end
