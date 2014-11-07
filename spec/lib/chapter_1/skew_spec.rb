require 'spec_helper'

xdescribe Bio::Skew do
  it 'calculates the skew' do
    skew = Bio::Skew.new('CATGGGCATCGGCCATACGCC')
    output = '0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2'
    expect(skew.find.join(' ')).to eq(output)
  end

  it 'finds the minimal skew positions' do
    skew = Bio::Skew.new('TAAAGACTGCCGAGAGGCCAACA' + 
                         'CGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT')
    expect(skew.minimal.join(' ')).to eq('11 24')
  end

end
