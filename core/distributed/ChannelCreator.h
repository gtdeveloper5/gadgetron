//
// Created by dchansen on 1/16/19.
//

#pragma once

#include <memory>
#include "../Channel.h"

namespace Gadgetron::Core::Distributed {
    class ChannelCreator {
    public:
        virtual std::shared_ptr<OutputChannel> create() = 0;
        virtual ~ChannelCreator() = default;
    };
}



