


int main(int argc, char *argv[]) {

  init_options();
  cmdline_getoptions();

  if (master_hostname) {
    fscl_slave(master_hostname, master_port);
    return 0;
  }

  verify_options();
  
  scan_obj = load_snp_input(snp_fname, include_invariant);
  
}
