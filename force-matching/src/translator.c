#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cgff_types.h"
#include "io_read.h"
#include "gromacs_topology.h"


int main(int argc, char * argv[])
{
  int ipt = -1, opt = -1, ft = -1, st = -1;
  char *ftc, *optc, *iptc, *stc;
  tW_word fnmi, fnmo, fnms, ns;
  bool bFT, bFI, bFO, bS, bIPT, bOPT, bST;
  strcpy(fnmi,"input.txt");
  strcpy(fnmo,"output.txt");
  strcpy(fnms,"N/A");
  strcpy(ns,"N/A");
  ftc = ns;
  optc = ns;
  iptc = ns;
  stc = ns;
  int i;

  tW_gmx_trxframe *fri = init_tW_gmx_trxframe();
  tW_gmx_trxframe *fro = init_tW_gmx_trxframe();
  tW_gmx_topology *top = init_tW_gmx_topology();



  bFT = bFI = bFO = bS = bIPT = bOPT = bST = FALSE ;


  void (*wtop)(tW_word, struct tW_gmx_topology *);

//parse args
  for (i = 0; i < argc; ++i)
  {
    if (strcmp(argv[i],"-ipt") == 0)
    {
      if (bIPT)
      {
	fprintf(stderr,"ERROR: found command line flag -ipt twice!\n");
        return 1;
      }
      bIPT = TRUE;
      check_count(i,argc,"-ipt");
      iptc = argv[i+1];
      if (strcmp(argv[i+1],DUMP) == 0) { ipt = fri->eFileType = top->eFileType = eDUMP;	}
      else if (strcmp(argv[i+1],BOCS) == 0) { ipt = fri->eFileType = top->eFileType = eBOCS; }
      else if (strcmp(argv[i+1],GRO) == 0) { ipt = fri->eFileType = eGRO; }
      else if (strcmp(argv[i+1],TRJ) == 0) { ipt = fri->eFileType = eTRJ; }
      else if (strcmp(argv[i+1],TRR) == 0) { ipt = fri->eFileType = eTRR; }
      else if (strcmp(argv[i+1],LMP) == 0) { ipt = fri->eFileType = eLMP; }
      else 
      {
	fprintf(stderr,"ERROR: Unrecognized value of -ipt: %s\n",argv[i+1]);
	fprintf(stderr,"Supported types: %s %s %s %s %s %s\n",DUMP,BOCS,GRO,TRJ,TRR,LMP);
	return 1;
      }
    }
    else if (strcmp(argv[i],"-opt") == 0)
    {
      if (bOPT)
      {
	fprintf(stderr,"ERROR: found command line flag -opt twice!\n");
	return 1;
      }
      bOPT = TRUE;
      check_count(i,argc,"-opt");
      optc = argv[i+1];
      if (strcmp(argv[i+1],DUMP) == 0) { opt = fro->eFileType = eDUMP; }
      else if (strcmp(argv[i+1],BOCS) == 0) { opt = fro->eFileType = eBOCS; wtop = write_bocs_top; }
      else if (strcmp(argv[i+1],GRO) == 0) { opt = fro->eFileType = eGRO; }
      else if (strcmp(argv[i+1],TRJ) == 0) { opt = fro->eFileType = eTRJ; }
      else if (strcmp(argv[i+1],TRR) == 0) { opt = fro->eFileType = eTRR; }
      else if (strcmp(argv[i+1],LMP) == 0) { opt = fro->eFileType = eLMP; }
      else
      {
	fprintf(stderr,"ERROR: Unrecognized value of -opt: %s\n",argv[i+1]);
	fprintf(stderr,"Supported types: %s %s %s %s %s %s\n",DUMP,BOCS,GRO,TRJ,TRR,LMP);
        return 1;
      }
    }
    else if (strcmp(argv[i],"-ft") == 0)
    {
      if (bFT)
      {
	fprintf(stderr,"ERROR: found command line flag -ft twice!\n");
	return 1;
      }
      bFT = TRUE;
      check_count(i,argc,"-ft");
      ftc = argv[i+1];
      if (strcmp(argv[i+1],TOP) == 0) { ft = eTOP; }
      else if (strcmp(argv[i+1],TRAJ) == 0) { ft = eTRAJ; }
      else
      {
	fprintf(stderr,"ERROR: Unrecognized value of -ft: %s\n",argv[i+1]);
	fprintf(stderr,"Supported types: %s %s\n",TOP,TRAJ);
	return 1;
      }
    }
    else if (strcmp(argv[i],"-st") == 0)
    {
      if (bST)
      {
        fprintf(stderr,"ERROR: found command line flag -st twice!\n");
        return 1;
      }
      bST = TRUE;
      check_count(i,argc,"-st");
      stc = argv[i+1];
      if (strcmp(argv[i+1],DUMP) == 0) { st = top->eFileType = eDUMP; }
      else if (strcmp(argv[i+1],BOCS) == 0) { st = top->eFileType = eBOCS; }
      else 
      { 
	fprintf(stderr,"ERROR: unrecognized option for -st: %s\n",argv[i+1]); 
	fprintf(stderr,"Supported types: %s %s\n",DUMP,BOCS); 
	return 1; 
      }
    }
    else if (strcmp(argv[i],"-fi") == 0) 
    { 
      if (bFI)
      {
        fprintf(stderr,"ERROR: found command line flag -fi twice!\n");
        return 1;
      }
      bFI = TRUE; 
      check_count(i,argc,"-fi"); 
      strcpy(fnmi,argv[i+1]); 
    }
    else if (strcmp(argv[i],"-fo") == 0) 
    {
      if (bFO)
      {
        fprintf(stderr,"ERROR: found command line flag -fo twice!\n");
        return 1;
      }
      bFO = TRUE; 
      check_count(i,argc,"-fo"); 
      strcpy(fnmo,argv[i+1]); 
    }  
    else if (strcmp(argv[i],"-s") == 0) 
    { 
      if (bS)
      {
        fprintf(stderr,"ERROR: found command line flag -s twice!\n");
        return 1;
      }
      bS = TRUE; 
      check_count(i,argc,"-s"); 
      strcpy(fnms,argv[i+1]); 
    }
  }
  
  if (ft == eTOP)
  {
    top->eFileType = ipt;
  }
  else
  {
    top->eFileType = st;
  }

  copyright();

// Print the values that were provided

  fprintf(stderr,"Option     Filename          Type         Description\n");
  fprintf(stderr,"--------------------------------------------------------------------\n");
  fprintf(stderr,"  -fi      %16s  Input        File to translate\n",fnmi);
  fprintf(stderr,"  -fo      %16s  Output       Translated file \n",fnmo);
  fprintf(stderr,"  -s       %16s  Input, Opt.  Topology file (req. if using -opt gro)\n",fnms);
  fprintf(stderr,"\n");
  fprintf(stderr,"Option       Type   Value   Description\n");
  fprintf(stderr,"------------------------------------------------------\n");
  fprintf(stderr,"-h           bool           Print help info and quit\n");
  fprintf(stderr,"-ft        string   %5s   Flag for what type of file you are translating (required)\n",ftc);
  fprintf(stderr,"                            options: %s %s\n",TRAJ,TOP);
  fprintf(stderr,"-ipt       string   %5s   Flag for format of -fi (required)\n",iptc);
  fprintf(stderr,"                            options: %s %s %s %s %s %s\n",DUMP,GRO,BOCS,TRJ,TRR,LMP);
  fprintf(stderr,"-opt       string   %5s   Flag for format of -fo (required)\n",optc);
  fprintf(stderr,"                            options: %s %s %s %s %s %s\n",DUMP,GRO,BOCS,TRJ,TRR,LMP);
  fprintf(stderr,"-st        string   %5s   Flag for format of -s (required if -s is specified)\n",stc);
  fprintf(stderr,"                            options: %s %s\n",DUMP,BOCS);

  for (i=0;i<argc;i++)
  {
    if (strcmp(argv[i],"-h") == 0)
    {
      return 0;
    }
  }

  // Check to make sure all the command line args are set

  if (!bFI)
  {
    fprintf(stderr,"ERROR: Please specify an input file with flag -fi\n");
    return 1;
  }
  if (!bFO)
  {
    fprintf(stderr,"ERROR: Please specify an output file with flag -fo\n");
    return 1;
  }
  if (!bFT)
  {
    fprintf(stderr,"ERROR: Please specify a file type with flag -ft\n");
    return 1;
  }
  if (!bIPT)
  {
    fprintf(stderr,"ERROR: Please specify an input file type with flag -ipt\n");
    return 1;
  }
  if (!bOPT)
  {
    fprintf(stderr,"ERROR: Please specify an output file type with flag -opt\n");
    return 1;
  }
  if (bS && !bST)
  {
    fprintf(stderr,"ERROR: You specified topology file: %s\nPlease specify a topology file type with flag -st\n",fnms);
    return 1;
  }


  if (ft == eTOP)
  {
    if (opt != eBOCS) { fprintf(stderr,"ERROR: topology can only be written in BOCS format\n"); return 1; }
    if ((ipt != eBOCS) && (ipt != eDUMP)) 
    { 
      fprintf(stderr,"ERROR: topology can only be read from DUMP or BOCS format\n"); 
      return 1; 
    }
    top->b_tpr = read_topology(top,fnmi);
    wtop(fnmo,top);
  }
  else if (ft == eTRAJ)
  {
    open_write_trajectory(fro,fnmo);
    read_first_frame(fri,fnmi);
    if ((opt == eGRO) || (opt == eLMP))
    {
      if (st == -1)
      {
	fprintf(stderr,"ERROR: -opt %s requires inclusion of a corresponding topology file!\n",((opt == eGRO) ? "gro" : "lmp"));
        fprintf(stderr,"Please specify a toplogy file with -s and a topology file type with -st\n");
        return 1;
      }
      top->b_tpr = read_topology(top,fnms);
      set_atom_labels(fri,top);
      if (top->contents->atoms.nr != fri->contents->natoms)
      {
	fprintf(stderr,"ERROR: topology file %s contains information on %d atoms\n",fnms,top->contents->atoms.nr);
	fprintf(stderr,"However, trajectory file %s contains information on %d atoms\n",fnmi,fri->contents->natoms);
	return 1;
      }
    }
    if (opt == eLMP)
    {
      read_next_frame(fri);
      fclose(fro->fp);
      write_lammps_data(fri,top,fnmo);
      free(fri);
      free(top);
      return 0;   
    }
    fro->contents = fri->contents;
    while (read_next_frame(fri))
    {
      copy_trxframe_info(fri,fro);  
      write_frame(fro);
    }   
    fprintf(stderr,"\n");
    fclose(fri->fp);
    fclose(fro->fp);
    //Program is about to end so this isn't necessary
    //BUT I spent a couple hours fixing a memory leak
    //in a different part of the program so right now
    //I'm free-memory crazy
    if (fri->contents->x) { free(fri->contents->x); }
    if (fri->contents->v) { free(fri->contents->v); }
    if (fri->contents->f) { free(fri->contents->f); }
    free(fri->contents);
    free(fri);
    free(fro);
    free(top->contents);
    free(top);
  }
  return 0;
}

