import React from 'react';
import PureRenderMixin from 'react-addons-pure-render-mixin';
import FieldGroupToStoreMixin from '../fieldGroup/FieldGroupToStoreMixin.js';

import InputFieldLabel from './InputFieldLabel.jsx';

var RadioGroupInputField= React.createClass(
    {
        mixins: [PureRenderMixin, FieldGroupToStoreMixin],

        propTypes: {
            options: React.PropTypes.array.isRequired,
        },


        componentWillMount() {
            // if no default value is specified, select the first option
            if (typeof this.state.fieldState.value === 'undefined' || this.state.fieldState.value==='') {
                this.state.fieldState.value = this.props.options[0].value;
            }
        },

        onChange(ev) {
            // the value of the group is the value of the selected option
            var val = ev.target.value;
            var checked = ev.target.checked;

            if (checked) {
                this.fireValueChange({
                    fieldKey: this.props.fieldKey,
                    newValue: val,
                    fieldState: this.state.fieldState
                });
            }
        },

        render() {
            return (
                <div style={{whiteSpace:'nowrap'}}>
                    <InputFieldLabel label={this.getLabel()}
                        tooltip={this.getTip()}
                        labelWidth={this.props.labelWidth}
                    />
                    {this.props.options.map((option) => {
                        var alignment = typeof option.align=='undefined'? 'top':option.align;//this is added by LZ
                        return (

                        <div key={option.value} style={{display:'inline-block'}}>
                            <input type='radio'
                                align={alignment }//this is added by LZ
                                name={this.props.fieldKey}
                                value={option.value}
                                defaultChecked={this.getValue()===option.value}
                                onChange={this.onChange}
                            />&nbsp;{option.label}&nbsp;&nbsp;
                        </div>

                        );
                        })}
                  </div>
            );
        }


    });

export default RadioGroupInputField;

