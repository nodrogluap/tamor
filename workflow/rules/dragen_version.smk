is_dragen_v42 = True
is_dragen_v44 = False
is_dragen_v45 = False
nirvana_downloader_path = "/opt/edico/share/nirvana/Downloader"
if(not os.path.exists(nirvana_downloader_path)):
        is_dragen_v42 = False
        # In Dragen v4.4+ it's changed to be from the PATH so that multiple Dragen versions can be supported on one server.
        result = subprocess.run(["/usr/bin/which", "dragen_info"], capture_output=True, text=True, check=True)
        if result.stderr:
                raise SystemExit("FATAL: Aborting, did not find Dragen install in path (assuming Dragen v4.4+)")
        install_path = result.stdout.removesuffix("/bin/dragen_info\n")
        if "4.5" in install_path:
                is_dragen_v45 = True
        elif "4.4" in install_path:
                is_dragen_v44 = True
        else:   
                raise SystemExit("FATAL: Unsupported Dragen version (not 4.2, 4.4 or 4.5) detected in dragen_info's install path (" + install_path+")")

