require 'spec_helper'

describe Bio::OriC do
  let(:vibrio)       { './spec/fixtures/vibrio-cholerae.txt' }
  let(:frequent)     { './spec/fixtures/frequent.txt' }
  let(:output_path)  { './spec/fixtures/output.txt' }
  let(:vibrio_dna)   { IO.readlines(vibrio).first }
  let(:frequent_dna) { IO.readlines(frequent).first }
  let(:output)       { IO.readlines(output_path).first }


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

end
