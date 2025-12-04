
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include "Datastructures/qlp/Qlp.hpp"
#include "yInterface.h"

void buildQlp(data::Qlp& qlp) {

    qlp.clear();

    /*// lia_1.qlp
    qlp.createVariable(data::QpVar("x", 0, 0.0, 3000.0,  data::QpVar::NumberSystem::generals, data::QpVar::all));
    qlp.createVariable(data::QpVar("y", 1, 0.0, 2001.0,  data::QpVar::NumberSystem::generals, data::QpVar::exists));    
    qlp.createVariable(data::QpVar("z", 2, 0.0, 1,  data::QpVar::NumberSystem::binaries, data::QpVar::exists));

    std::vector<data::IndexedElement> r1;
    r1.emplace_back(0, 2);
    r1.emplace_back(1, -1);
    r1.emplace_back(2, -10000);
    data::Constraint& c1e = qlp.createRhsConstraint(data::QpRhs(-1, data::QpRhs::smallerThanOrEqual,data::Constraint::EXISTENTIAL));
    c1e.setElements(r1);
    std::vector<data::IndexedElement> r2;
    r2.emplace_back(0, -1);
    r2.emplace_back(1, -1);
    r2.emplace_back(2, +10000);
    data::Constraint& c2e = qlp.createRhsConstraint(data::QpRhs(6998,data::QpRhs::smallerThanOrEqual,data::Constraint::EXISTENTIAL));
    c2e.setElements(r2);  
*/

    std::vector<data::QpVar> x = qlp.createVariables("x", 4, data::QpVar::NumberSystem::binaries, 0.0, 1, data::QpVar::exists);
    data::QpVar y = qlp.createVariable("y", data::QpVar::NumberSystem::generals, 0.0, 10.0, data::QpVar::exists, -10);    
    data::QpVar y_bar =qlp.createVariable("y_bar", data::QpVar::NumberSystem::generals, 0.0, 10.0, data::QpVar::all);           
   

    qlp.setObjective(data::QpObjFunc::min);
    qlp.setObjectiveFunctionElement(x[0], -1);
    qlp.setObjectiveFunctionElement(x[1], -2);
    qlp.setObjectiveFunctionElement(x[2], -4);
    qlp.setObjectiveFunctionElement(x[3], -8);

    // -25 * x + 20 * y <= 30
    data::Constraint& c1e = qlp.createRhsConstraint(data::QpRhs::smallerThanOrEqual, 30, data::Constraint::EXISTENTIAL);
    for (int i=0;i<x.size();i++)
      qlp.addConstraintElement(c1e, x[i], -25*pow(2,i));
    qlp.addConstraintElement(c1e, "y", 20);

    data::Constraint& c1u = qlp.createRhsConstraint(data::QpRhs::smallerThanOrEqual, 30, data::Constraint::UNIVERSAL);
    for (int i=0;i<x.size();i++)
       qlp.addConstraintElement(c1u, x[i], -25*pow(2,i));
    qlp.addConstraintElement(c1u, y_bar, 20);
    
    // x + 2 * y <= 10
    data::Constraint& c2e = qlp.createRhsConstraint(data::QpRhs::smallerThanOrEqual, 10, data::Constraint::EXISTENTIAL);
    for (int i=0;i<x.size();i++)
      qlp.addConstraintElement(c2e, x[i], 1*pow(2,i));
    qlp.addConstraintElement(c2e, y, 2);

    data::Constraint& c2u = qlp.createRhsConstraint(data::QpRhs::smallerThanOrEqual, 10, data::Constraint::UNIVERSAL);
    for (int i=0;i<x.size();i++)
      qlp.addConstraintElement(c2u, x[i], 1*pow(2,i));
    qlp.addConstraintElement(c2u, y_bar, 2);

   
    // 2 * x - y <= 15
    data::Constraint& c3e = qlp.createRhsConstraint(data::QpRhs::smallerThanOrEqual, 15, data::Constraint::EXISTENTIAL);
    for (int i=0;i<x.size();i++)
      qlp.addConstraintElement(c3e, x[i], 2*pow(2,i));
    qlp.addConstraintElement(c3e, y, -1);
    data::Constraint& c3u = qlp.createRhsConstraint(data::QpRhs::smallerThanOrEqual, 15, data::Constraint::UNIVERSAL);
    for (int i=0;i<x.size();i++)
      qlp.addConstraintElement(c3u, x[i], 2*pow(2,i));
    qlp.addConstraintElement(c3u, y_bar, -1);

    
    // 2 * x + 10 * y >= 15
    data::Constraint& c4e = qlp.createRhsConstraint(data::QpRhs::greaterThanOrEqual, 15, data::Constraint::EXISTENTIAL);
    for (int i=0;i<x.size();i++)
      qlp.addConstraintElement(c4e, x[i], 2*pow(2,i));
    qlp.addConstraintElement(c4e, y, 10);
    data::Constraint& c4u = qlp.createRhsConstraint(data::QpRhs::greaterThanOrEqual,15, data::Constraint::UNIVERSAL);
    for (int i=0;i<x.size();i++)
      qlp.addConstraintElement(c4u, x[i], 2*pow(2,i));
    qlp.addConstraintElement(c4u, y_bar, 10);

    //Binarize x
    data::Constraint& binX = qlp.createRhsConstraint(data::QpRhs::smallerThanOrEqual,10, data::Constraint::EXISTENTIAL);
    for (int i=0;i<x.size();i++)
      qlp.addConstraintElement(binX, x[i], pow(2,i));

    // y - y_bar <= 0
    data::Constraint& c5e = qlp.createRhsConstraint(data::QpRhs::smallerThanOrEqual, 0, data::Constraint::EXISTENTIAL);
    qlp.addConstraintElement(c5e, y, 1);
    qlp.addConstraintElement(c5e, y_bar, -1);

    std::cout << qlp.toString() << std::endl;
}

int main(int argc, char** argv) {
  data::Qlp* orgQlpPt = new data::Qlp;
  buildQlp(*orgQlpPt);  // Build your instance using the QLP datastructure
  
  yInterface* yasolPt=new yInterface(); // Instantiate Yasol
  yasolPt->yInit(*orgQlpPt);            // Initialize Yasol
  //yasolPt->setParam("showWarning",true); 
  //yasolPt->setParam("useGMI",3); 
  //yasolPt->setTimelimit(20);             // Set Timelimit (in seconds)
  //yasolPt->setGap(20);             // Set Gap (in percent)

  yasolPt->supressOutput();
  yasolPt->solve();                     // Call Yasol 
  //yasolPt->writeSolutionFile("Out.sol");
  
  std::cout << "Solution Status: " << yasolPt->getStatus() << std::endl;
  std::cout << "RESULT: " << yasolPt->getResult() << std::endl;
  std::cout << "GAP: " << yasolPt->getGap() << std::endl;
  std::cout << "DUAL: " << yasolPt->getDual() << std::endl;

  std::cout << "TIME: " << yasolPt->getRuntime() << std::endl;
  if(yasolPt->getStatus()==YASOL_OPTIMAL || yasolPt->getStatus()==YASOL_FEASIBLE || yasolPt->getStatus()==YASOL_INCUMBENT){
    vector<SolutionEntry> theSolution= yasolPt->getSolution();
    int block=0;
    for (int i=0;i<theSolution.size();i++){
      if (block<theSolution[i].block){
        block++;
        cout << "==============Block "<< block << "=============="<<endl; 
      }
      cout << theSolution[i].name<<" = "<< theSolution[i].value << endl;
    }
  }

  delete orgQlpPt;
  delete yasolPt;
  return 0;
}

