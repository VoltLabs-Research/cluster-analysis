import json
import os

from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps, cmake_layout


class ClusterAnalysisConan(ConanFile):
    name = "cluster-analysis"
    package_type = "static-library"
    license = "MIT"
    settings = "os", "arch", "compiler", "build_type"
    default_options = {"hwloc/*:shared": True}
    requires = (
        "boost/1.88.0",
        "onetbb/2021.12.0",
        "coretoolkit/[>=2.5]",
        "spdlog/1.14.1",
        "nlohmann_json/3.11.3",
    )
    exports = "vpm.json"
    exports_sources = "CMakeLists.txt", "include/*", "src/*"

    def set_version(self):
        manifest = os.path.join(self.recipe_folder, "vpm.json")
        with open(manifest, encoding="utf-8") as stream:
            self.version = json.load(stream)["version"]

    def layout(self):
        cmake_layout(self)

    def generate(self):
        CMakeToolchain(self).generate()
        CMakeDeps(self).generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.set_property(
            "cmake_target_name", "cluster-analysis::cluster-analysis"
        )
        self.cpp_info.libs = ["cluster-analysis_lib"]
        self.cpp_info.requires = [
            "boost::headers",
            "onetbb::onetbb",
            "coretoolkit::coretoolkit",
            "nlohmann_json::nlohmann_json",
            "spdlog::spdlog",
        ]
