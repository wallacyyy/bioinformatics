require 'spec_helper'

describe Bio::PatternCounter do
  let(:sample)  { 'GCGCG' }
  let(:pattern) { 'GCG' }
  let(:dna) do
    File.open('./spec/fixtures/vibrio-cholerae.txt', 'rb') { |file| file.read }
  end

  it 'should count overlapping patterns' do
    counter = Bio::PatternCounter.new(sample, pattern)
    expect(counter.count).to eq(2)
  end

  it 'should count patterns from a normal size dataset in a reasanable time' do
    expect(Benchmark.realtime{ Bio::PatternCounter.new(dna, 'TGA').count } < 1).to eq(true)
  end

end
