#ifndef KinematicLV_hh
#define KinematicLV_hh
class KineamticLV: public TLorentzVector{
	private:
		double Mass = 0;
		vector<TLorentzVector> D;
	public:
		KinematicLV(double M_, vector<TLorentzVector> D_);
}
#endif
