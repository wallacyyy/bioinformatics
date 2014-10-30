require 'spec_helper'

describe Bio::OriC do
  let(:vibrio_dna)   { IO.read('./spec/fixtures/vibrio-cholerae.txt') }
  let(:frequent_dna) { IO.read('./spec/fixtures/frequent.txt') }
  let(:output)       { IO.read('./spec/fixtures/output.txt') }

  it 'counts overlapping patterns' do
    sample  = 'GCGCG'
    pattern = 'GCG'

    counter = Bio::OriC.new(sample)
    expect(counter.count(pattern)).to eq(2)
  end

  it 'counts patterns from a normal size dataset in a reasonable time' do
    benchmark = Benchmark.realtime{ Bio::OriC.new(vibrio_dna).count('TGA') }
    expect(benchmark < 1).to eq(true)
  end

  context 'finds the most common pattern given a k-mer' do
    it 'runs with a simple dna' do
      sample = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
      k_mer  = 4
      result = Bio::OriC.new(sample).most_frequents(k_mer)

      expect(result).to eq(['CATG', 'GCAT'])
    end

    it 'runs with a complex dna' do
      k_mer = 12
      result = Bio::OriC.new(frequent_dna).most_frequents(k_mer)
      expect(result).to eq(output.split)
    end
  end

  it 'finds frequency tables' do
    sample = IO.read('./spec/fixtures/frequency-table.txt').chomp
    output = IO.read('./spec/fixtures/frequency-table-output.txt').chomp
    bio = Bio::OriC.new(sample)
    expect(bio.frequency_table(6).join(' ')).to eq(output)
  end

  context 'find clumps' do
    it 'runs with a simple dna string' do
      bio = Bio::OriC.new('CGGACTCGACAGATGTGAAGAACGACAATGTG' +
                          'AAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA')
      expect(bio.clumps(5, 50, 4)).to eq(['CGACA', 'GAAGA'])
    end

    xit 'runs with a great dna string' do
      sample = IO.read('./spec/fixtures/clumps.txt').chomp
      bio = Bio::OriC.new(sample)
      expect(bio.clumps(11, 566, 18)).to eq(['AAACCAGGTGG'])
    end

  end
end
