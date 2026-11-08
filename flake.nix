{
  description = "MUSIC pole-zero mismatch figures";
  inputs = {
    nixpkgs.follows = "utils/nixpkgs";
    flake-utils.url = "github:numtide/flake-utils";
    utils.url = "github:ewtodd/Analysis-Utilities";
    # The MUSIC monorepo, pinned. Only its tooling PulseHistory, Constants and
    # Plots are compiled into the binaries; nothing there is edited.
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
        analysis-utils = utils.packages.${system}.default;
        root = utils.packages.${system}.root;
        musicTooling = "${music.outPath}/tooling";

        # Keep simulation output and build artifacts out of the store copy.
        src = lib.cleanSourceWith {
          src = ./.;
          filter =
            path: _type:
            let
              rel = lib.removePrefix (toString ./. + "/") (toString path);
              top = lib.head (lib.splitString "/" rel);
            in
            !(lib.elem top [
              "output"
              "doc"
              "bin"
              "build"
            ])
            && !lib.hasPrefix "result" top;
        };
      in
      {
        # `nix build` -> ./result/bin/figures
        packages.default = pkgs.stdenv.mkDerivation {
          pname = "polezero-music-figures";
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

        devShells.default = pkgs.mkShell {
          nativeBuildInputs = with pkgs; [
            pkg-config
            gnumake
            clang-tools
          ];
          buildInputs = [
            analysis-utils
            root
            pkgs.bash
          ];
          shellHook = ''
            echo "Analysis-Utilities version: ${analysis-utils.version}"
            export MUSIC_TOOLING="''${MUSIC_TOOLING:-${musicTooling}}"
            echo "MUSIC tooling: $MUSIC_TOOLING"
            export CPLUS_INCLUDE_PATH="$PWD/include:$MUSIC_TOOLING/include''${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"
            export ROOT_INCLUDE_PATH="$PWD/include:$MUSIC_TOOLING/include''${ROOT_INCLUDE_PATH:+:$ROOT_INCLUDE_PATH}"
            alias clean-aclic='rm -f *_C.so *_C.d *_C_ACLiC_dict_rdict.pcm *_cpp.so *_cpp.d *_cpp_ACLiC_dict_rdict.pcm AutoDict_*'
          '';
        };
      }
    );
}
