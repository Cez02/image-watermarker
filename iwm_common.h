#pragma once

#include <string>
#include <stdexcept>
#include <format>

namespace iwm {
    
    enum class IWM_Mode {
        lsb_spatial,
        pvd_spatial,
        ae_spatial,
        dft_me,
        dct_fe,
        midb_emb,
        MODE_COUNT
    };

    struct IWM_ModeDescription {
        std::string Code, FullName;
    };

    const IWM_ModeDescription IWM_MODEDESCRIPTIONS[static_cast<uint>(IWM_Mode::MODE_COUNT)] = {
        { "lsb_spatial", "Least Significant Bit Insertion" },
        { "pvd_spatial", "Pixel Value Differencing" },
        { "ae_spatial", "Additive Embedding" },
        { "dft_me", "DFT Magnitude Embedding (LSB)" },
        { "dct_fe", "DCT Frequency Embedding (LSB)" },
        { "midb_emb", "Middle band DFT-DCT Embedding" },
    };

    IWM_Mode parse_iwmmode(std::string str){
        uint i;

        for(i = 0; i<static_cast<uint>(IWM_Mode::MODE_COUNT); i++){
            if(IWM_MODEDESCRIPTIONS[i].Code == str){
                return static_cast<IWM_Mode>(i);
            }
        }

        throw std::invalid_argument(std::format("{} is not a recognizable IWM Mode.", str));
    }

}