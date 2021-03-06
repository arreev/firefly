/*
 * License information at https://github.com/Caltech-IPAC/firefly/blob/master/License.txt
 */
import validator from 'validator';
import {isNil} from 'lodash';

var isInRange= function(val,min,max) {
    var retval= !(min !== undefined && min!==null && val<min);
    return retval && !(max !== undefined && max!==null && val>max);
};

var typeInject= {
    asInt : {
        dataTypeDesc : 'integer',
        testFunc : validator.isInt,
        convertFunc : validator.toInt
    },
    asFloat : {
        dataTypeDesc : 'float',
        testFunc : validator.isFloat,
        convertFunc : validator.toFloat
    }
};



var makePrecisionStr= function(value,precision) {
    if (!isNil(value)) {
        return (precision>-1) ? value.toFixed(precision) : value;
    }
    else return '';
};

var makeErrorMessage= function(description,min,max,precision) {
    var retval= '';
    var hasMin= (min !== undefined && min!==null);
    var hasMax= (min !== undefined && min!==null);
    var minStr= makePrecisionStr(min,precision);
    var maxStr= makePrecisionStr(max,precision);
    description= description || '';
    if (hasMin && hasMax) {
        retval= description + ': must be between ' + minStr + ' and ' + maxStr;
    }
    else if (hasMin) {
        retval= description + ': must be greater than ' + minStr;
    }
    else if (hasMax) {
        retval= description + ': must be less than ' + maxStr;
    }
    return retval;
};

var validateRange = function(min,max,precision,description,dType, valStr, nullAllowed) {
    var retval= {
        valid : true,
        message : ''
    };
    if (valStr) {
        valStr+= '';
        if (valStr && dType.testFunc(valStr)) {
            const v = dType.convertFunc(valStr);
            if (!isInRange(v, min, max) || isNaN(v)) {
                retval.valid = false;
                retval.message = makeErrorMessage(description, min, max, precision);
            }
        }
        else {
            retval.valid = false;
            retval.message = description + ': must be a '+ dType.dataTypeDesc + makeErrorMessage(null, min, max, precision);
        }
    }
    else if (!nullAllowed) {
        retval.valid = false;
        retval.message = description + ': must be a '+ dType.dataTypeDesc + makeErrorMessage(null, min, max, precision);
    }
    return retval;
};

export const validateEmail = function(description,valStr) {
    var retval = {
        valid: true,
        message: ''
    };
    if (valStr && !validator.isEmail(valStr+'')) {
        retval.valid = false;
        retval.message = description + ': must be a valid email address';
    }
    return retval;
};

export const validateUrl = function(description,valStr) {
    var retval = {
        valid: true,
        message: ''
    };
    if (valStr && !validator.isURL(valStr+'')) {
        retval.valid = false;
        retval.message = description + ': must be a valid URL';
    }
    return retval;
};

export const intRange = function(min,max,description, valStr, nullAllowed=true) {
   return validateRange(min,max,null,description,typeInject.asInt,valStr, nullAllowed);
};

export const floatRange = function(min,max,precision, description, valStr, nullAllowed=true) {
    return validateRange(min,max,precision,description,typeInject.asFloat,valStr, nullAllowed);
};

export const isFloat = function(description, valStr) {
    var retval= { valid : true, message : '' };
    if (valStr) {
        if (!validator.isFloat(valStr+'')) {
            retval.valid = false;
            retval.message = description + ': must be a float';
        }
    }
    return retval;
};

export const isPositiveFiniteNumber = (description, valStr)=>{
    var retval= { valid : true, message : '' };
    if (valStr) {
        var aNumber = Number.parseFloat(valStr);
        if (!isFinite(aNumber)) {
            retval.valid = false;
            retval.message = description + ': must be a finite float';
        }
        if (aNumber<0){
            retval.valid = false;
            retval.message = description + ': must be a positive float';
        }
    }
    return retval;

};


export const isInt = function(description, valStr) {
    var retval= { valid : true, message : '' };
    if (valStr) {
        if (!validator.isInt(valStr+'')) {
            retval.valid = false;
            retval.message = description + ': must be an int';
        }
    }
    return retval;
};

export const isHexColorStr = function(description, valStr) {
    var retval= { valid : true, message : '' };
    if (valStr && !/^#[0-9a-f]{6}/.test(valStr)) {
        retval.valid = false;
        retval.message = description + ': must be a hex color exactly 7 characters long';
    }

};


/*---------------------------- validator function used by InputField to validate a value -----------------------------*/
/*---- these factory functions creates a validation function that takes a value and return {valid,message} -----------*/
export const intValidator = function(min,max,description) {
    return (val) => intRange(min, max, description, val);
};

export const floatValidator = function(min,max,description) {
    return (val) => floatRange(min, max, description, val);
};

export const urlValidator = function(description) {
    return (val) => validateUrl(description, val);
};

export const emailValidator = function(description) {
    return (val) => validateEmail(description, val);
};
/*--------------------------------------------------------------------------------------------------------------------*/



var Validate = {
    validateEmail, validateUrl, intRange, floatRange, isFloat, isInt,isPositiveFiniteNumber
};
export default Validate;