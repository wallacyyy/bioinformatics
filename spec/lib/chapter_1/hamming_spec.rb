require 'spec_helper'

xdescribe Bio::Hamming do
  let(:hamming) { Bio::Hamming.new }

  it 'find frequent patterns by sort' do
    sample = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
    result = hamming.frequent_sort(sample, 4)
    expect(result).to eq(['CATG', 'GCAT'])
  end

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
      result = hamming.frequent_patterns_sort(sample, 4, 1)
      expect(result.join(' ')).to eq('ATGC ATGT GATG')
    end

    it 'runs with a complex dna string' do
      sample = IO.read('./spec/fixtures/mismatch.txt').chomp
      result = hamming.frequent_patterns_sort(sample, 10, 2)
      expect(result.join(' ')).to eq('GCACACAGAC GCGCACACAC')
    end

    it 'finds it considering their reverse complements' do
      sample = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
      result = hamming.frequent_patterns_sort(sample, 4, 1, true)
      expect(result.join(' ')).to eq('ACAT ATGT')
    end

    it 'finds it considering their reverse complements 
        with a complex dna string' do
      sample = IO.read('./spec/fixtures/mismatch-reverse.txt').chomp
      result = hamming.frequent_patterns_sort(sample, 9, 3, true)
      expect(result.join(' ')).to eq('AGCGCCGCT AGCGGCGCT')
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
