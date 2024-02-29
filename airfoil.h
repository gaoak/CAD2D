#ifndef AIRFOIL_H
#define AIRFOIL_H
#include <map>
#include <memory>
#include <string>
#include <vector>
class AirFoil;

typedef std::shared_ptr<AirFoil> AirFoilShPtr;

class AirFoil {
public:
  AirFoil() {}
  virtual std::vector<double> Upper(double x) = 0;
  virtual std::vector<double> Lower(double x) = 0;
  virtual void GetInfo(std::map<std::string, double> &p) = 0;
};

class NACAmpxx : public AirFoil {
public:
  NACAmpxx(std::map<std::string, double> params);
  NACAmpxx(double m, double p, double t, bool roundtrailing = false);
  NACAmpxx(std::string name);
  std::vector<double> Upper(double x) override;
  std::vector<double> Lower(double x) override;
  double Findx(double s, int up = 1);
  double Finds(double x, int up = 1);
  void GetInfo(std::map<std::string, double> &p) override;

protected:
  void InitAirfoil();
  double Findx(double s, std::vector<std::vector<double>> &arc);
  double Finds(double x, std::vector<std::vector<double>> &arc);
  double halfSt(double x);
  double halft(double x);
  double halfdt(double x);
  double chamber(double x);
  double theta(double x);
  double calculateTrailingRadius(double &xtmp);
  double testRadius(double x, double &r);
  void calculateArcTable();
  std::vector<double> roundTrailingEdge(std::vector<double> &p0,
                                        double eps = 1.E-8);
  double roundTrailingSize();
  double Area();
  std::vector<std::vector<double>> m_arcu;
  std::vector<std::vector<double>> m_arcd;
  double m_m;
  double m_p;
  double m_t;
  double m_TERadius;
  double m_TETangencyX;
  bool m_isRoundTrailing;
};

class WedgeFoil : public AirFoil {
public:
  WedgeFoil(std::map<std::string, double> params);
  WedgeFoil(double D0, double D1, double D2, bool roundtrailing = false);
  std::vector<double> Upper(double x) override;
  std::vector<double> Lower(double x) override;
  double Findx(double s, int up = 1);
  double Finds(double x, int up = 1);
  void GetInfo(std::map<std::string, double> &p) override;

protected:
  void InitAirfoil();
  std::vector<double> roundTrailingEdge(std::vector<double> &p0,
                                        double eps = 1.E-8);
  double roundTrailingSize();
  double Area();
  double m_LEDiameter;
  double m_TEThich;
  double m_TERadius;
  double m_theta; // polar angle of the tangency point betweent the leading-edge
                  // circule and the wedge, < 90 degrees
  double
      m_LETangencyX; // 0.5 * m_LEDiameter + 0.5 * m_LEDiameter * cos(m_theta)
  double m_TETangencyX; // 1. - m_TEDiameter * (1 - cos(m_theta))
  bool m_isRoundTrailing;
};

#endif
