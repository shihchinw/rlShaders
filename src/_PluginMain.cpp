#include <string>
#include <ai.h>

extern  AtNodeMethods *GgxMethod;
extern  AtNodeMethods *DisneyMethod;

enum    ShaderId
{
    kGgx,
    kDisney,
};


node_loader
{
    switch (i) {
        case kGgx:
            node->methods = GgxMethod;
            node->output_type = AI_TYPE_RGB;
            node->name = "rlGgx";
            node->node_type = AI_NODE_SHADER;
            break;

        case kDisney:
            node->methods = DisneyMethod;
            node->output_type = AI_TYPE_RGB;
            node->name = "rlDisney";
            node->node_type = AI_NODE_SHADER;
            break;

        default:
            return false;
    }
    
    strcpy_s(node->version, AI_VERSION);

    return true;
}