// $Id: inline_pion_decay.cc, v1.2 2015-07-05 22:29:00 sdcohen Exp $
/*! \file
 * \brief Inline construction of pion-decay two-point correlators
 *
 * Calculates the two-point correlator of the pseudoscalar meson corresponding
 * to leptonic decay by temporal-axial current.
 *
 * v1.2 adds input of momentum using a list
 * v1.1 adds multiple-source location trickery
 */

//#include "inline_eigen_maker.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "util/info/proginfo.h"
#include "io/param_io.h"
#include "io/qprop_io.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "overlap_eigenop.h"

namespace Chroma {

//------------------------------------------------------------------------------
// Factory Registration

// Environment in which the measurement lives (holds params and such)
namespace InlineEigenMakerEnv {
namespace {
// Function to register with a factory
AbsInlineMeasurement *createMeasurement(XMLReader &xml_in,
                                        const std::string &path) {
  return new InlineEigenMaker(InlineEigenMakerParams(xml_in, path));
}
AbsInlineMeasurement *createMeasurement2(XMLReader &xml_in,
                                         const std::string &path) {
  return new InlineEigenEraser(InlineEigenEraserParams(xml_in, path));
}

//! Local registration flag
bool registered = false;
}

// The name of the measurement for the XML file
const std::string name = "EIGEN_MAKER";
const std::string name2 = "EIGEN_ERASER";

//! Register all the factories
bool registerAll() {
  bool success = true;
  if (!registered) {
    success &= TheInlineMeasurementFactory::Instance().registerObject(
        name, createMeasurement);
    success &= TheInlineMeasurementFactory::Instance().registerObject(
        name2, createMeasurement2);
    registered = true;
  }
  return success;
}
} // end namespace InlineEigenMakerEnv

//------------------------------------------------------------------------------
// Parameter reading, writing, handling

int check_eigen_name(std::string &name, int size) {
  // TODO: check whether the file exist with the corect size;
  FILE *pFile;
  long fsize;
  pFile = fopen(name.c_str(), "rb");
  if (pFile == NULL)
    return false;
  else {
    fseek(pFile, 0, SEEK_END);
    fsize = ftell(pFile);
    fclose(pFile);
    return true;
  }
}

//! Reader for parameters
void read(XMLReader &xml, const std::string &path,
          InlineEigenMakerParams::Param_t &param) {
  XMLReader paramtop(xml, path);

  param.fermact = readXMLGroup(paramtop, "FermionAction", "FermAct");
  read(paramtop, "Noeigen", param.Noeigen);
  read(paramtop, "gpu_level", param.gpu_level);
  read(paramtop, "check_residual", param.check_residual);

  param.filename = "_NULL_";
  if (paramtop.count("filename") > 0)
    read(paramtop, "filename", param.filename);

  QDPIO::cout << param.filename << std::endl;

  param.file_exist = false;
  // check whether the file exist with reasonable size
  if (param.filename != "_NULL_")
    param.file_exist = check_eigen_name(param.filename, param.Noeigen);

  if (param.file_exist == false) {
    read(paramtop, "extra_space", param.extra_space);
    read(paramtop, "chebyshev_cut", param.chebyshev_cut);
    read(paramtop, "chebyshev_order", param.chebyshev_order);
    read(paramtop, "use_ckpoint", param.use_ckpoint);
  } else {
    param.extra_space = 0;
    param.chebyshev_cut = 0.3;
    param.chebyshev_order = 0;
    param.use_ckpoint = false;
  }
}

//! Writer for parameters
void write(XMLWriter &xml, const std::string &path,
           const InlineEigenMakerParams::Param_t &param) {
  push(xml, path);

  push(xml, "FermionAction");
  xml << param.fermact.xml;
  pop(xml);
  write(xml, "Noeigen", param.Noeigen);
  write(xml, "gpu_level", param.gpu_level);
  write(xml, "check_residual", param.check_residual);

  write(xml, "filename", param.filename);

  write(xml, "extra_space", param.extra_space);
  write(xml, "chebyshev_cut", param.chebyshev_cut);
  write(xml, "chebyshev_order", param.chebyshev_order);
  write(xml, "use_ckpoint", param.use_ckpoint);

  pop(xml);
}

//! Named objects (gauge config, prop pairs) input
void read(XMLReader &xml, const std::string &path,
          InlineEigenMakerParams::NamedObject_t &input) {
  XMLReader inputtop(xml, path);

  read(inputtop, "gauge_id", input.gauge_id);
  read(inputtop, "eigen_id", input.eigen_id);
  read(inputtop, "op_id", input.op_id);
}

//! Named objects (gauge config, prop pairs) output
void write(XMLWriter &xml, const std::string &path,
           const InlineEigenMakerParams::NamedObject_t &input) {
  push(xml, path);

  write(xml, "gauge_id", input.gauge_id);
  write(xml, "eigen_id", input.eigen_id);
  write(xml, "op_id", input.op_id);

  pop(xml);
}

// Params default constructor
InlineEigenMakerParams::InlineEigenMakerParams(void) { frequency = 0; }

// Construct params from XML
InlineEigenMakerParams::InlineEigenMakerParams(XMLReader &xml_in,
                                               const std::string &path) {
  try {
    XMLReader paramtop(xml_in, path);

    if (paramtop.count("Frequency") == 1)
      read(paramtop, "Frequency", frequency);
    else
      frequency = 1;

    // Read in the parameters
    read(paramtop, "Param", param);

    // Read in the output propagator/source configuration info
    read(paramtop, "NamedObject", named_obj);

    // Possible alternate XML file pattern
    if (paramtop.count("xml_file") != 0) {
      read(paramtop, "xml_file", xml_file);
    }
  }
  catch (const std::string &e) {
    QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e
                << std::endl;
    QDP_abort(1);
  }
}

// Write out the parameters we constructed
void InlineEigenMakerParams::write(XMLWriter &xml_out,
                                   const std::string &path) {
  push(xml_out, path);

  Chroma::write(xml_out, "Param", param);
  Chroma::write(xml_out, "NamedObject", named_obj);
  QDP::write(xml_out, "xml_file", xml_file);

  pop(xml_out);
}

//------------------------------------------------------------------------------
// Helper functions to read in pairs of propagators

// Anonymous namespace for propagator-handling utilities
//------------------------------------------------------------------------------
// The real work is done here

// Set up the XML and invoke func, which does the actual work
void InlineEigenMaker::operator()(unsigned long update_no, XMLWriter &xml_out) {
  // If xml file not empty, then use alternate
  if (params.xml_file != "") {
    std::string xml_file = makeXMLFileName(params.xml_file, update_no);

    push(xml_out, "EigenMaker");
    write(xml_out, "update_no", update_no);
    write(xml_out, "xml_file", xml_file);
    pop(xml_out);

    XMLFileWriter xml(xml_file);
    func(update_no, xml);
  } else {
    func(update_no, xml_out);
  }
}

// Real work done here
void InlineEigenMaker::func(unsigned long update_no, XMLWriter &xml_out) {
  START_CODE();

  StopWatch snoop;
  snoop.reset();
  snoop.start();

  // Boilerplate stuff to the output XML
  push(xml_out, "EigenMaker");
  write(xml_out, "update_no", update_no);

  // Write some diagnostic junk about the lattice
  QDPIO::cout << " EigenMaker: Make eigen system of " << params.param.fermact.id
              << " action." << std::endl;
  QDPIO::cout << "     volume: " << Layout::lattSize()[0];
  for (int i = 1; i < Nd; ++i) {
    QDPIO::cout << " x " << Layout::lattSize()[i];
  }
  QDPIO::cout << std::endl;

  // Write info about the program
  proginfo(xml_out);

  // Write out the input
  params.write(xml_out, "Input");

  EigenOperator<LatticeFermion> *eigen;
  OverlapEigenOperator *eigen_ov;

  push(xml_out, "EigenMaker_measurements");

  // if(TheNamedObjMap::Instance().check(params.named_obj.eigen_id))
  if (TheNamedObjMap::Instance().check(params.named_obj.op_id)) {
    QDPIO::cerr
        << params.named_obj.op_id
        << " has been assigned to some other object and can not be used here"
        << std::endl;
    QDP_abort(1);
  } else {
    if (params.param.fermact.id != "HWILSON" &&
        params.param.fermact.id != "OVERLAP") {
      QDPIO::cerr << "Action " << params.param.fermact.id
                  << " is not supported yet" << std::endl;
      QDP_abort(1);
    }

    TheNamedObjMap::Instance().create<EigenOperator<LatticeFermion> *>(
        params.named_obj.eigen_id);
    TheNamedObjMap::Instance().create<OverlapEigenOperator *>(
        params.named_obj.op_id);

    // Write out the config info
    // Grab gauge configuration
    multi1d<LatticeColorMatrix> U;
    XMLBufferWriter gauge_xml;
    try {
      // Try to get the gauge field.
      // If it doesn't exist, an exception will be thrown.
      U = TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix> >(
          params.named_obj.gauge_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(
          gauge_xml);
    }
    catch (std::bad_cast) {
      // If the object exists and is of the wrong type we will land here
      QDPIO::cerr << InlineEigenMakerEnv::name << ": caught dynamic cast error"
                  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string &e) {
      // If the object is not in the map we will land here
      QDPIO::cerr << InlineEigenMakerEnv::name << ": map call failed: " << e
                  << std::endl;
      QDP_abort(1);
    }
    // Bind gauge field reference to a local name
    const multi1d<LatticeColorMatrix> &u =
        TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix> >(
            params.named_obj.gauge_id);
    write(xml_out, "Config_info", gauge_xml);

    if (params.param.fermact.id == "HWILSON")
      eigen = new HwilsonEigenOperator(params.param.fermact, U,
                                       params.param.Noeigen);
    if (params.param.fermact.id == "OVERLAP") {
      eigen = new OverlapEigenOperator(params.param.fermact, U,
                                       params.param.Noeigen);

      TheNamedObjMap::Instance().getData<OverlapEigenOperator *>(
          params.named_obj.op_id) = eigen;
    }
    TheNamedObjMap::Instance().getData<EigenOperator<LatticeFermion> *>(
        params.named_obj.eigen_id) = eigen;

    // TheNamedObjMap::Instance().getData<OverlapEigenOperator* >
    //            (params.named_obj.op_id)=eigen;

    // Keep an array of all the xml output buffers
    if (params.param.file_exist == true) {
      eigen->load(params.param.filename);
    } else {
      // TODO:create the eigensystem;
      eigen->create_eigen(params);
      QDPIO::cout << "to save.." << std::endl;
      exit(0);
      std::string name = params.param.filename + ".test";
      eigen->save(name);
    }
  }

  if (params.param.check_residual == true) {
    eigen->check_residual(fuzz);
    QDPIO::cout << "Check for overlap operator" << std::endl;
  }

  pop(xml_out); // EigenMaker

  snoop.stop();
  QDPIO::cout << InlineEigenMakerEnv::name
              << ": total time = " << snoop.getTimeInSeconds() << " secs"
              << std::endl;

  QDPIO::cout << InlineEigenMakerEnv::name << ": ran successfully" << std::endl;

  END_CODE();
} // end of InlineEigenMaker::func

}; // end of namespace Chroma
