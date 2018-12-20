
#include <memory>

#include "connection/ProtoConnection.h"
#include "connection/ConfigConnection.h"
#include "connection/StreamConnection.h"
#include "connection/Writers.h"

#include "Connection.h"


#include "readers/Primitives.h"
#include "log.h"

using namespace boost::asio;

using namespace Gadgetron::Core;
using namespace Gadgetron::Core::Readers;

using namespace Gadgetron::Server::Connection;
using namespace Gadgetron::Server::Connection::Writers;

namespace {

    class ErrorChannel : public ErrorHandler {
    public:
        ErrorChannel() : errors(std::make_shared<MessageChannel>()) {}

        void handle(const std::string &location, std::function<void()> fn) override {
            try {
                fn();
            }
            catch (const std::exception &e) {
                push_error(location, e.what());
            }
            catch (const std::string &s) {
                push_error(location, s);
            }
            catch (...) {
                push_error(location, "Unknown exception.");
            }
        }

    private:
        void push_error(const std::string location, const std::string &message) {
            errors->push(std::make_unique<std::string>("[" + location + "] ERROR: " + message));
        }

        std::shared_ptr<MessageChannel> errors;
    };


    void send_errors(std::iostream &stream) {}

    void send_close(std::iostream &stream) {
        uint16_t close = 4;
        stream.write(reinterpret_cast<char *>(&close), sizeof(close));
    }

    void handle_connection(const Gadgetron::Core::Context::Paths &paths, std::unique_ptr<std::iostream> stream) {

        ErrorChannel error_handler{};

        auto config = ProtoConnection::process(*stream, paths, error_handler);

        if (config) {
            auto context = ConfigConnection::process(*stream, paths, error_handler);

            StreamConnection::process(*stream, context, config.get(), error_handler);
        }

        send_errors(*stream);

        send_close(*stream);
    }
}

namespace Gadgetron::Server::Connection {

    void handle(const Gadgetron::Core::Context::Paths &paths, std::unique_ptr<std::iostream> stream) {
        auto thread = std::thread(handle_connection, paths, std::move(stream));
        thread.detach();
    }

    Connection::Connection(std::iostream &stream) : stream(stream) {};

    void Connection::process_input() {

        bool closed = false;
        auto handlers = prepare_handlers(closed);

        while (!closed) {
            auto id = read_t<uint16_t>(stream);

            GDEBUG_STREAM("Processing message with id: " << id);

            handlers.at(id)->handle(stream);
        }
    }

    void Connection::process_output() {

        auto writers = prepare_writers();

        writers.push_back(std::make_unique<TextWriter>());
        writers.push_back(std::make_unique<ResponseWriter>());

        InputChannel<Message> &output = *channels.output;
        for (auto message : output) {

            auto writer = std::find_if(writers.begin(), writers.end(),
                   [&](auto &writer) { return writer->accepts(*message); }
            );

            (*writer)->write(stream, std::move(message));
        }
    }

    std::vector<std::unique_ptr<Writer>> Connection::prepare_writers() {
        return std::vector<std::unique_ptr<Writer>>();
    }

    void Connection::start(ErrorHandler &handler) {

        std::function<void()> handled_input = [&]() {
            handler.handle("Connection Input", [&]() { this->process_input(); });
        };

        std::function<void()> handled_output = [&]() {
            handler.handle("Connection Output", [&]() { this->process_output(); });
        };

        threads.input  = std::thread(handled_input);
        threads.output = std::thread(handled_output);
    }

    void Connection::join() {
        threads.input.join();
        threads.output.join();
    }
}
