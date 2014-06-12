// This file is released into the public domain

if getos() == "Windows" then
    // to manage long pathname
    includes_src_c = "-I""" + get_absolute_file_path("builder_gateway_c.sce") + "../../src/c""";
else
    includes_src_c = "-I" + get_absolute_file_path("builder_gateway_c.sce") + "../../src/c";
end

tbx_build_gateway("serialize_c", ..
                  ["serialize_set","sci_serialize_set";"serialize_get","sci_serialize_get"], ..
                  ["sci_serialize.c"], ..
                  get_absolute_file_path("builder_gateway_c.sce"), ..
                  ["../../src/c/libserialize"], ..
                  "", ..
                  includes_src_c);

clear tbx_build_gateway;
