{
  description = "MUSIC pole-zero mismatch simulation (L/R preamp decay-time study)";
  inputs = {
    nixpkgs.follows = "utils/nixpkgs";
    flake-utils.url = "github:numtide/flake-utils";
    utils = {
      url = "github:ewtodd/Analysis-Utilities";
    };
    # The MUSIC monorepo, pinned. Nothing of its own is built: the unmodified
    # tooling/src/PulseHistory.cpp and Constants.cpp are compiled into these
    # binaries from the input's source tree. Following utils' nixpkgs keeps the
    # whole tree on one pin; MUSIC's CUDA ROOT is never evaluated into these
    # derivations, so its binary-cache caveat does not apply here.
    music = {
      url = "github:ewtodd/MUSIC";
      inputs.nixpkgs.follows = "utils/nixpkgs";
    };
  };
  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      utils,
      music,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };
        lib = pkgs.lib;
        # CPU build of the library and the matching ROOT it was built against.
        # The simulation has no GPU work, so the CUDA variants are not needed.
        analysis-utils = utils.packages.${system}.default;
        root = pkgs.root;

        # $(MUSIC_TOOLING) of the Makefile: include/ and src/ under the pinned
        # input. Already a store path, so the build reads it directly; no copy,
        # no second src.
        musicTooling = "${music.outPath}/tooling";

        # Keep the tree out of the store copy: output/ holds ~60 MB of root
        # files, doc/ only matters to the scripts, and bin/, build/ and the
        # result symlink are build artifacts of `nix develop` and `nix build`.
        src = lib.cleanSourceWith {
          src = ./.;
          filter = path: _type:
            let
              rel = lib.removePrefix (toString ./. + "/") (toString path);
              top = lib.head (lib.splitString "/" rel);
            in
              !(lib.elem top [
                "output"
                "doc"
                "bin"
                "build"
              ]) && !lib.hasPrefix "result" top;
        };
      in
      {
        # `nix build` -> ./result/bin/{simulate,analyze,realwave,realdt,
        # kernelstudy,taustudy,figures,wavecheck,deconvolve,optfilter,
        # applicability}. Nothing is written into the tree.
        packages.default = pkgs.stdenv.mkDerivation {
          pname = "polezero-music-sim";
          version = "1.0";
          inherit src;

          nativeBuildInputs = with pkgs; [
            pkg-config
            gnumake
          ];
          buildInputs = [
            analysis-utils
            root
            pkgs.bash
          ];

          # The Makefile's MUSIC_TOOLING ?= falls back to ~/MUSIC/tooling;
          # this env value wins.
          MUSIC_TOOLING = toString musicTooling;

          dontConfigure = true;
          buildPhase = ''
            runHook preBuild
            make -j$NIX_BUILD_CORES
            runHook postBuild
          '';
          installPhase = ''
            runHook preInstall
            mkdir -p $out/bin
            cp bin/* $out/bin/
            chmod +x $out/bin/*
            runHook postInstall
          '';
        };

        # Interactive shell (ROOT/CLING work, editing against a local MUSIC
        # checkout with MUSIC_TOOLING=~/MUSIC/tooling). The canonical build is
        # `nix build`; the scripts in scripts/ use it too.
        devShells.default = pkgs.mkShell {
          nativeBuildInputs = with pkgs; [
            pkg-config
            gnumake
            clang-tools
            # Firmware reverse-engineering of the CAEN V2740 DPP-PHA .CUP
            # (ARM64 control-plane decompilation + device-tree clock inspection).
            ghidra
            dtc
          ];
          buildInputs = [
            analysis-utils
            root
            pkgs.bash
          ];
          shellHook = ''
            echo "Analysis-Utilities version: ${analysis-utils.version}"
            flake_root="$PWD"
            # Unmodified MUSIC tooling sources (PulseHistory and the Constants
            # it depends on). The pinned flake input is the default so a bare
            # `make` matches `nix build`; point MUSIC_TOOLING at a local
            # checkout to compile against working changes.
            export MUSIC_TOOLING="''${MUSIC_TOOLING:-${musicTooling}}"
            echo "MUSIC tooling: $MUSIC_TOOLING"
            export CPLUS_INCLUDE_PATH="$flake_root/include:$MUSIC_TOOLING/include''${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"
            export ROOT_INCLUDE_PATH="$flake_root/include:$MUSIC_TOOLING/include''${ROOT_INCLUDE_PATH:+:$ROOT_INCLUDE_PATH}"
            alias clean-aclic='rm -f *_C.so *_C.d *_C_ACLiC_dict_rdict.pcm *_cpp.so *_cpp.d *_cpp_ACLiC_dict_rdict.pcm AutoDict_*'
          '';
        };
      }
    );
}
