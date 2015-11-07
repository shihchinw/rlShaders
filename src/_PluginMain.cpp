#include <string>
#include <ai.h>

extern  AtNodeMethods *GgxMethod;
extern  AtNodeMethods *DisneyMethod;
//extern  AtNodeMethods *DisneyXMethod;
extern  AtNodeMethods *SssNDMethod;
extern  AtNodeMethods *SkinMethod;

enum    ShaderId
{
    kGgx,
    kDisney,
    //kDisneyX,
    kSssND,
    kSkin,
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

       /* case kDisneyX:
            node->methods = DisneyXMethod;
            node->output_type = AI_TYPE_RGB;
            node->name = "rlDisneyX";
            node->node_type = AI_NODE_SHADER;
            break;*/

        case kSssND:
            node->methods = SssNDMethod;
            node->output_type = AI_TYPE_RGB;
            node->name = "rlSssND";
            node->node_type = AI_NODE_SHADER;
            break;

        case kSkin:
            node->methods = SkinMethod;
            node->output_type = AI_TYPE_RGB;
            node->name = "rlSkin";
            node->node_type = AI_NODE_SHADER;
            break;

        default:
            return false;
    }

    strcpy_s(node->version, AI_VERSION);

    return true;
}