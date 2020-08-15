/*! \file
 * \brief Inline construction of propagator
 *
 * PropagatorMulti calculations
 */

#include "fermact.h"
#include "inline_propagator_multi_eigen.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "util/info/unique_id.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlinePropagatorMultiEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlinePropagatorMulti(InlinePropagatorMultiParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "PROPAGATOR_MULTI";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= WilsonTypeFermActsEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  } // end namespace


  //! PropagatorMulti input
  void read(XMLReader& xml, const std::string& path, InlinePropagatorMultiParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "source_id", input.source_id);
    read(inputtop, "prop_id", input.prop_id);
  }

  //! PropagatorMulti output
  void write(XMLWriter& xml, const std::string& path, const InlinePropagatorMultiParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "source_id", input.source_id);
    write(xml, "prop_id", input.prop_id);

    pop(xml);
  }


  // Param stuff
  InlinePropagatorMultiParams::InlinePropagatorMultiParams() { frequency = 0; }

  InlinePropagatorMultiParams::InlinePropagatorMultiParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Parameters for source construction
      read(paramtop, "Param", param);

      read(paramtop, "mass", mass);

      // Read in the output propagator/source configuration info
      read(paramtop, "NamedObject", named_obj);
      
      if(named_obj.prop_id.size()!=mass.size())
      {
         QDPIO::cerr << __func__ << "Size of propagators doesn't equal to that of the masses" << std::endl;
         QDP_abort(1);
      }

      // Possible alternate XML file pattern
      if (paramtop.count("xml_file") != 0) 
      {
	read(paramtop, "xml_file", xml_file);
      }
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
      QDP_abort(1);
    }
  }


  void
  InlinePropagatorMultiParams::writeXML(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    write(xml_out, "Param", param);
    write(xml_out, "Param", mass);
    write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }


  // Function call
  void 
  InlinePropagatorMulti::operator()(unsigned long update_no,
			       XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      std::string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "propagator");
      write(xml_out, "update_no", update_no);
      write(xml_out, "xml_file", xml_file);
      pop(xml_out);

      XMLFileWriter xml(xml_file);
      func(update_no, xml);
    }
    else
    {
      func(update_no, xml_out);
    }
  }

  void quarkPropMult(multi1d<LatticePropagator>& psi,const multi1d<Real>& masses,
                        const LatticePropagator& chi, 
                        int j_decay,
                        const ChromaProp_t& Param_,
                        int &ncg_had)
  {
  
  
  }

  // Real work done here
  void 
  InlinePropagatorMulti::func(unsigned long update_no,
			 XMLWriter& xml_out) 
  {
    START_CODE();

    QDPIO::cout << InlinePropagatorMultiEnv::name << ": propagator calculation" << std::endl;

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    // Test and grab a reference to the gauge field
    XMLBufferWriter gauge_xml;
    try
    {
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlinePropagatorMultiEnv::name << ": caught dynamic cast error" 
		  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << InlinePropagatorMultiEnv::name << ": std::map call failed: " << e 
		  << std::endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "propagator");
    write(xml_out, "update_no", update_no);

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.writeXML(xml_out, "Input");

    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    //
    // Read in the source along with relevant information.
    // 
    XMLReader source_file_xml, source_record_xml;

    // These pesky variables are needed in the quarkprop call - only chiral dudes
    // need this stuff, but it must be there for the bleeping virtual function
    // to live at the base class
    int t0;
    int j_decay;

    // Record the type of header
    bool make_sourceP = false;
    bool seqsourceP = false;

    QDPIO::cout << "Snarf the source from a named buffer" << std::endl;
    try
    {
      // Try the cast to see if this is a valid source
      LatticePropagator& source_tmp = 
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.source_id);

      // Snarf the source info. This is will throw if the source_id is not there
      TheNamedObjMap::Instance().get(params.named_obj.source_id).getFileXML(source_file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.source_id).getRecordXML(source_record_xml);

      // Try to invert this record XML into a source struct
      // First identify what kind of source might be here
      if (source_record_xml.count("/MakeSource") != 0)
      {
	make_sourceP = true;
	MakeSourceProp_t  orig_header;
	read(source_record_xml, "/MakeSource", orig_header);

	j_decay = orig_header.source_header.j_decay;
	t0      = orig_header.source_header.t_source;
      }
      else if (source_record_xml.count("/SequentialSource") != 0)
      {
	seqsourceP = true;
	SequentialSource_t   orig_header;
	read(source_record_xml, "/SequentialSource", orig_header);

	j_decay = orig_header.seqsource_header.j_decay;
	t0      = orig_header.seqsource_header.t_sink;   // funky, but probably not really needed
      }
      else
      {
	throw std::string("No appropriate header found");
      }

      // Write out the source header
      write(xml_out, "Source_file_info", source_file_xml);
      write(xml_out, "Source_record_info", source_record_xml);
    }    
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlinePropagatorMultiEnv::name << ": caught dynamic cast error" 
		  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << InlinePropagatorMultiEnv::name << ": error extracting source_header: " << e << std::endl;
      QDP_abort(1);
    }

    // Should be a valid cast now
    const LatticePropagator& quark_prop_source = 
      TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.source_id);
 
    QDPIO::cout << "Source successfully read and parsed" << std::endl;

    // Sanity check - write out the norm2 of the source in the Nd-1 direction
    // Use this for any possible verification
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, Nd-1);

      multi1d<Double> source_corr = sumMulti(localNorm2(quark_prop_source), 
					     phases.getSet());

      push(xml_out, "Source_correlator");
      write(xml_out, "source_corr", source_corr);
      pop(xml_out);
    }

    //
    // Loop over the source color and spin, creating the source
    // and calling the relevant propagator routines. The QDP
    // terminology is that a propagator is a matrix in color
    // and spin space
    //
    multi1d<LatticePropagator > quark_propagator(0);
    try
    {
      quark_propagator.resize(params.mass.size());
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlinePropagatorMultiEnv::name << ": caught dynamic cast error" 
		  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << InlinePropagatorMultiEnv::name << ": error creating prop: " << e << std::endl;
      QDP_abort(1);
    }

    int ncg_had = 0;

    //
    // Initialize fermion action
    //
    std::istringstream  xml_s(params.param.fermact.xml);
    XMLReader  fermacttop(xml_s);
    QDPIO::cout << "FermAct = " << params.param.fermact.id << std::endl;

    //
    // Try the factories
    //
    bool success = false;

    if (! success)
    {
      try
      {
	StopWatch swatch;
	swatch.reset();
	QDPIO::cout << "Try the various factories" << std::endl;

	   QDPIO::cout << "Compute the multiple quark props" << std::endl;
	   swatch.start();

           quarkPropMult(quark_propagator, 
                       params.mass, 
                       quark_prop_source,
                       j_decay,
                       params.param, 
                       ncg_had);
	
           swatch.stop();
	   QDPIO::cout << "PropagatorMulti computed: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << std::endl;

	success = true;
      }
      catch (const std::string& e) 
      {
	QDPIO::cout << InlinePropagatorMultiEnv::name << ": caught exception around quarkprop: " << e << std::endl;
      }
    }


    if (! success)
    {
      QDPIO::cerr << "Error: no fermact found" << std::endl;
      QDP_abort(1);
    }


    push(xml_out,"Relaxation_Iterations");
    write(xml_out, "ncg_had", ncg_had);
    pop(xml_out);

    // Sanity check - write out the propagator (pion) correlator in the Nd-1 direction
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, Nd-1);

      for(int im=0;im<params.mass.size();im++)
      {
         multi1d<Double> prop_corr = sumMulti(localNorm2(quark_propagator[im]), 
					   phases.getSet());

         push(xml_out, "Prop_correlator");
         write(xml_out, "prop_corr", prop_corr);
         pop(xml_out);
      }
    }


    // Save the propagator info
    try
    {
      QDPIO::cout << "Start writing propagator info" << std::endl;

      XMLBufferWriter file_xml;
      push(file_xml, "propagator");
      write(file_xml, "id", uniqueId());  // NOTE: new ID form
      pop(file_xml);

      XMLBufferWriter record_xml;
      if (make_sourceP)
      {
	MakeSourceProp_t  orig_header;
	read(source_record_xml, "/MakeSource", orig_header);

	Propagator_t new_header;   // note, abandoning state_info
	new_header.prop_header   = params.param;
	new_header.source_header = orig_header.source_header;
	new_header.gauge_header  = orig_header.gauge_header;
	write(record_xml, "PropagatorMulti", new_header);  
      } 
      else if (seqsourceP)
      {
	SequentialSource_t  orig_header;
	read(source_record_xml, "/SequentialSource", orig_header);

	SequentialProp_t  new_header;   // note, abandoning state_info
	new_header.seqprop_header   = params.param;
	new_header.sink_header      = orig_header.sink_header;
	new_header.seqsource_header = orig_header.seqsource_header;
	new_header.forward_props    = orig_header.forward_props;
	new_header.gauge_header     = orig_header.gauge_header;
	write(record_xml, "SequentialProp", new_header);  
      }

      // Write the propagator xml info
      for(int im=0;im<params.mass.size();im++)
      {
          TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.prop_id[im]);
          TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_id[im])=quark_propagator[im];
          TheNamedObjMap::Instance().get(params.named_obj.prop_id[im]).setFileXML(file_xml);
          TheNamedObjMap::Instance().get(params.named_obj.prop_id[im]).setRecordXML(record_xml);
      }

      QDPIO::cout << "PropagatorMulti successfully updated" << std::endl;
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlinePropagatorMultiEnv::name << ": caught dynamic cast error" 
		  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << InlinePropagatorMultiEnv::name << ": error extracting prop_header: " << e << std::endl;
      QDP_abort(1);
    }

    pop(xml_out);  // propagator

    snoop.stop();
    QDPIO::cout << InlinePropagatorMultiEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << std::endl;

    QDPIO::cout << InlinePropagatorMultiEnv::name << ": ran successfully" << std::endl;

    END_CODE();
  } 

}
